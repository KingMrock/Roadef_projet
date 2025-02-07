import json
from docplex.mp.model import Model
import argparse
import time
import os

def compute_tight_M(instance):
    """
    Compute a tight big-M for the risk constraints.
    For each intervention i and time t, we loop over every valid starting time s (that makes i active at t)
    and take the maximum risk (over the risk scenarios). Then for each time t, we sum the maximum risk
    contributions over interventions and finally take the maximum over t.
    """
    interventions = instance["Interventions"]
    T = int(instance["T"])
    # R[i][t] holds the maximum risk for intervention i at time t.
    R = { i: { t: 0 for t in range(1, T+1) } for i in interventions }


    
    for i, data in interventions.items():
        tmax = int(data["tmax"])
        Delta = data["Delta"]  # list of durations for each possible start time
        risk_data = data.get("risk", {})
        for t in range(1, T+1):
            for s in range(1, tmax+1):
                # intervention i is active at time t when started at s if:
                if s <= t < s + Delta[s-1]:
                    # Get risk list for time t and start s (as a list over scenarios)
                    risk_list = risk_data.get(str(t), {}).get(str(s), [])
                    if risk_list:
                        R[i][t] = max(R[i][t], max(risk_list))
    # For each time period, sum the maximum risk contributions from all interventions.
    M_ts = { t: sum(R[i][t] for i in interventions) for t in range(1, T+1) }
    M_risk = max(M_ts.values())
    return M_risk

def risk_value(i, t, s, interventions, n_scenarios):
    """
    Return the risk list for intervention i at time t if it starts at s.
    If the risk list is missing or too short, it is padded with zeros.
    """
    risk_data = interventions[i].get("risk", {})
    risk_list = risk_data.get(str(t), {}).get(str(s), [])
    if len(risk_list) < n_scenarios:
        risk_list = risk_list + [0]*(n_scenarios - len(risk_list))
    return risk_list


def build_model(instance):
    mdl = Model("Grid_Interventions")
    
    interventions = instance["Interventions"]
    T = int(instance["T"])
    resources = instance.get("Resources", {})
    seasons = instance.get("Seasons", {})
    
    # Objective parameters
    alpha = float(instance.get("Alpha", 0.5))
    quantile = float(instance.get("Quantile", 0.5))
    
    # Compute a tight M for the risk constraints.
    M_risk = compute_tight_M(instance)
    print("Computed tight M for risk constraints:", M_risk)

    print("Building Decision Variables")
    
    # 1. Decision Variables: For each intervention i and each possible start s (1..tmax)
    x = {}
    for i, data in interventions.items():
        tmax = int(data["tmax"])
        x[i] = {}
        for s in range(1, tmax+1):
            x[i][s] = mdl.binary_var(name=f"x_{i}_{s}")
        mdl.add_constraint(
            mdl.sum(x[i][s] for s in range(1, tmax+1)) == 1,
            ctname=f"start_once_{i}"
        )

    print("Adding Ressource Constraints")
    
    # 2. Resource Constraints
    for r, r_data in resources.items():
        max_array = r_data.get("max", [999999]*T)
        min_array = r_data.get("min", [-999999]*T)
        resource_usage = { t: [] for t in range(1, T+1) }
        for i, data in interventions.items():
            tmax_i = int(data["tmax"])
            wdict_res = data.get("workload", {}).get(r, {})
            for s in range(1, tmax_i+1):
                for t_str, s_dict in wdict_res.items():
                    t_int = int(t_str)
                    usage = s_dict.get(str(s), 0)
                    if usage != 0:
                        resource_usage[t_int].append((i, s, usage))
        for t in range(1, T+1):
            if resource_usage[t]:
                usage_expr = mdl.sum(x[i][s] * usage for (i, s, usage) in resource_usage[t])
                mdl.add_constraint(usage_expr <= max_array[t-1], ctname=f"capacity_{r}_{t}")
                mdl.add_constraint(usage_expr >= min_array[t-1], ctname=f"min_capacity_{r}_{t}")
    

    print("Adding Active Variables and Exclusions")
    # 3. Active Variables for Exclusions: active[i,t]=1 if intervention i is active at time t.
    active = {}
    for i, data in interventions.items():
        active[i] = {}
        tmax_i = int(data["tmax"])
        delta_list = data["Delta"]
        for t in range(1, T+1):
            a_var = mdl.binary_var(name=f"active_{i}_{t}")
            active[i][t] = a_var
            cover_expr = []
            for s in range(1, tmax_i+1):
                dur_s = delta_list[s-1]
                if s <= t <= s + dur_s - 1:
                    cover_expr.append(x[i][s])
            mdl.add_constraint(a_var == mdl.sum(cover_expr), ctname=f"def_active_{i}_{t}")
    
    # 3b. Exclusion Constraints: For each exclusion [i1, i2, season], at times in that season.
    for ex_name, ex_list in instance.get("Exclusions", {}).items():
        i1, i2, season_name = ex_list
        if season_name in seasons:
            for t_str in seasons[season_name]:
                t_int = int(t_str)
                mdl.add_constraint(
                    active[i1][t_int] + active[i2][t_int] <= 1,
                    ctname=f"excl_{i1}_{i2}_{t_int}"
                )
    

    print("Compute Risk Model")
    # 4. Risk Model
    # Determine number of scenarios (assume risk lists have the same length).
    scenario_count = None
    for i, data in interventions.items():
        if "risk" in data:
            for t_str, s_dict in data["risk"].items():
                for s_str, risk_list in s_dict.items():
                    scenario_count = len(risk_list)
                    break
                if scenario_count is not None:
                    break
        if scenario_count is not None:
            break
    if scenario_count is None:
        scenario_count = 1
    scenario_indices = range(scenario_count)
    

    print("Computing Quantile values")
    # Create Q[t] variables and risk expression (and z[t,sc] binary flags).
    Q = { t: mdl.continuous_var(name=f"Q_{t}") for t in range(1, T+1) }
    z = {}
    risk_expr = {}
    for t in range(1, T+1):
        for sc in scenario_indices:
            # Sum contributions from all interventions that are active at time t.
            expr = mdl.sum(
                x[i][s] * risk_value(i, t, s, interventions, scenario_count)[sc]
                for i, data in interventions.items()
                for s in range(1, int(data["tmax"])+1)
                if s <= t < s + data["Delta"][s-1]
            )
            risk_expr[(t, sc)] = expr
            z[(t, sc)] = mdl.binary_var(name=f"z_{t}_{sc}")
            mdl.add_constraint(
                expr <= Q[t] + M_risk * z[(t, sc)],
                ctname=f"risk_bound_{t}_{sc}"
            )
        mdl.add_constraint(
            mdl.sum(z[(t, sc)] for sc in scenario_indices) <= (1 - quantile) * scenario_count,
            ctname=f"quantile_{t}"
        )
    

    print("Building Objective Function")
    # Compute mean risk over scenarios at each time t.
    mean_risk = {}
    for t in range(1, T+1):
        risk_expressions = [risk_expr[(t, sc)] for sc in scenario_indices]
        mean_risk[t] = (1.0 / scenario_count) * mdl.sum(risk_expressions)
    
    # 5. Excess variables: e[t] >= Q[t] - mean_risk[t] and e[t] >= 0.
    e = {}
    for t in range(1, T+1):
        e[t] = mdl.continuous_var(lb=0, name=f"excess_{t}")
        mdl.add_constraint(e[t] >= Q[t] - mean_risk[t], ctname=f"def_excess_{t}")
    
    # 6. Objective: weighted combination of mean risk and excess.
    obj1_expr = (1.0 / T) * mdl.sum(mean_risk[t] for t in range(1, T+1))
    obj2_expr = (1.0 / T) * mdl.sum(e[t] for t in range(1, T+1))
    mdl.minimize(alpha * obj1_expr + (1 - alpha) * obj2_expr)

    
    # Store x in the model so we can retrieve the decision values later.
    mdl._x = x
    return mdl

    

def main():
    
    parser = argparse.ArgumentParser(description="Solve the intervention scheduling problem.")
    parser.add_argument("instance_file", type=str, help="Path to the JSON instance file.")
    parser.add_argument("--solution", "-s", type=str, default="example_solution.txt", help="Path to the solution export file.")
    parser.add_argument("--time", "-t", type=int, default=600, help="Time limit for the solver in seconds.")
    
    args = parser.parse_args()

    if not os.path.exists(args.instance_file):
        print(f"Instance file {args.instance_file} not found.")
        return
    
    with open(args.instance_file, "r") as f:
        instance = json.load(f)
    
    interventions = instance["Interventions"]
    
    #get start time
    start = time.time()

    print("Building model...")

    mdl = build_model(instance)


    #get time needed to build the model
    curr = time.time()-start

    print("Model built in ", curr, " seconds")

    #set time limit for the solver as the remaining time
    args.time = args.time - curr

    #print remaining time
    print("Time remaining for the solver: ", args.time, " seconds")


    if args.time:
        mdl.parameters.timelimit = args.time
    
    # Optionally, adjust solver parameters, e.g., MIP gap tolerance
    mdl.parameters.mip.tolerances.mipgap = 0.01  # 1% gap tolerance

    # Print the number of variables and constraints
    print("Number of variables:", mdl.number_of_variables)
    print("Number of constraints:", mdl.number_of_constraints)
    
    solution = mdl.solve(log_output=True)
    if solution:
        print("Objective value:", solution.objective_value)
        with open(args.solution, "w") as fsol:
            for i, data_i in interventions.items():
                tmax_i = int(data_i["tmax"])
                start_chosen = None
                for s in range(1, tmax_i+1):
                    if solution.get_value(mdl._x[i][s]) > 0.5:
                        start_chosen = s
                        break
                fsol.write(f"{i} {start_chosen}\n")
    else:
        print("No solution found.") 


if __name__ == "__main__":
    main()
