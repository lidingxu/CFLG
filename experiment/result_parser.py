import os
import copy
from typing import NamedTuple
import collections
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt

min_primal_bound = 2**31
max_dual_bound = -1
time_limit = 3600
algorithms  = ["EF", "EFP", "LEFP", "EVFP", "LEVFP", "None"] # ["LEFPI", "LEFP", "LEFPD", "LEFPV",  "None"] # ["LEFPI", "LEFP", "LEFPD", "LEFPV", "None"]
coverages = ["Small", "Large"]
benchmarks = ["city", "Kgroup_A", "Kgroup_B", "random_A", "random_B"]


comp_algos = ["LEFPD", "LEFPI", "LEFP", "LEFPV"]

styles = ['solid', 'dashed', 'dotted', 'dashdot']
markers = ['o', "d", "*", "X" ]
colors = ['magenta', 'blue', 'red', 'green']

algonm_map = {"EF": "EF", "LEFPI" : "EF-PLI", "LEFP" : "EF-PLB", "LEFPD": "EF-PLD", "LEFPV": "EF-PLBC", "LEFPV2": "EF-PLBC2",  "EVFP" : "EVF-PB", "LEVFP" : "EVF-PLB", "EFP":"EF-PB", "EFPD": "EF-PD", "EFPV": "EF-PV1", "EFPV2": "EF-PV2"}

def parse_name(name):
    s = ""
    for c in name:
        if c == "_":
            s += "\_"
        else:
            s +=c
    return s

def extractInstanceResult(file_path):
    ls = open(file_path).readlines()
    #print(file)
    stat_keys = ["time", "node", "gap", "obj", "bound", "instance", "formulation"]
    stat_dict = {}
    entries = {}
    for l in ls:
        for stat_key in stat_keys:
            if  l.split(":")[0].strip() == stat_key:
                stat_dict[stat_key] = l
    #print(stat_dict["soltime"].split(),  stat_dict["node"].split(), stat_dict["obj"].split(), stat_dict["relgap"].split())
    entries["time"] = float(stat_dict["time"].split()[1])
    entries["node"] = int(stat_dict["node"].split()[1])
    entries["obj"] = float(stat_dict["obj"].split()[1])
    entries["bound"] = float(stat_dict["bound"].split()[1])
    entries["gap"] = abs(entries["obj"] - entries["bound"])/ entries["obj"] * 100# float(stat_dict["relgap"].split()[1])*10000
    #print(entries["gap"])
    entries["absgap"] = abs(entries["obj"] - entries["bound"])
    entries["instance"] = stat_dict["instance"].split()[1]
    entries["formulation"] = stat_dict["formulation"].split()[1]
    entries["missing"] = False
    if entries["absgap"] > 1e20 or entries["bound"] > 1e20:
        entries["missing"] = True
        entries["gap"] = 100
    #print(entries, file_path)

    return entries
    #writer.writerow([name, algo, primal_bound, gap, total_time, pricing_time, pricing_calls, columns_gen, nodes])

def extractInstance(file_path):
    ls = open(file_path).readlines()
    #print(file)
    stat_keys = ["instance", "formulation", "dlt", "org_node", "org_edge", "org_min_len", "org_avg_len", "org_max_len",  "sdb_node", "sdb_edge", "sdb_min_len", "sdb_avg_len", "sdb_max_len",  "dtf_node", "dtf_edge", "dtf_min_len", "dtf_avg_len", "dtf_max_len"]
    stat_dict = {}
    entries = {}
    for l in ls:
        #print(l)
        for stat_key in stat_keys:
            if  l.split(":")[0].strip() == stat_key:
                stat_dict[stat_key] = l
    #print(stat_dict)
    #print(stat_dict["soltime"].split(),  stat_dict["node"].split(), stat_dict["obj"].split(), stat_dict["relgap"].split())
    entries["dlt"] = 1 #float(stat_dict["dlt"].split()[1])
    entries["org_node"] = float(stat_dict["org_node"].split()[1]) / entries["dlt"]
    entries["org_edge"] = float(stat_dict["org_edge"].split()[1]) / entries["dlt"]
    entries["org_min_len"] = float(stat_dict["org_min_len"].split()[1]) / entries["dlt"]
    entries["org_avg_len"] = float(stat_dict["org_avg_len"].split()[1]) / entries["dlt"]
    entries["org_max_len"] = float(stat_dict["org_max_len"].split()[1]) / entries["dlt"]
    entries["sdb_node"] = float(stat_dict["sdb_node"].split()[1]) / entries["dlt"]
    entries["sdb_edge"] = float(stat_dict["sdb_edge"].split()[1]) / entries["dlt"]
    entries["sdb_min_len"] = float(stat_dict["sdb_min_len"].split()[1]) / entries["dlt"]
    entries["sdb_avg_len"] = float(stat_dict["sdb_avg_len"].split()[1]) / entries["dlt"]
    entries["sdb_max_len"] = float(stat_dict["sdb_max_len"].split()[1]) / entries["dlt"]
    entries["dtf_node"] = float(stat_dict["dtf_node"].split()[1]) / entries["dlt"]
    entries["dtf_edge"] = float(stat_dict["dtf_edge"].split()[1]) / entries["dlt"]
    entries["dtf_min_len"] = float(stat_dict["dtf_min_len"].split()[1]) / entries["dlt"]
    entries["dtf_avg_len"] = float(stat_dict["dtf_avg_len"].split()[1]) / entries["dlt"]
    entries["dtf_max_len"] = float(stat_dict["dtf_max_len"].split()[1]) / entries["dlt"]
    entries["instance"] = stat_dict["instance"].split()[1]
    entries["formulation"] = stat_dict["formulation"].split()[1]
    entries["isnotfind"] = False
    #print(entries)
    return entries


defualt_entry = {"gap": 100.0, "absgap": 100.0, "time":1800, "node": 0,  "cdual": 100.0, "cprimal": 100.0, "obj": 0,  "obj_": 0, "type": 0}

shift = { "gap": 1.0, "time": 1.0, "solved": 0, "node": 1.0, "cdual": 1.0, "cprimal": 1.0, "obj": 1, "obj_": 1}
sgm_keys = ["gap", "time", "node", "cdual", "cprimal", "obj", "obj_"]
display_keys = [ "time", "solved",  "total", "node", "gap", "obj"]
percent_keys = ["gap", "obj", "obj_"]
int_keys = [ "solved", "total"]


def Stat(algo, coverage):
    return {"algorithm": algo, "cover": coverage, "solved": 0, "total": 0, "solution": 0, "gap": 0.0, "time":0.0, "node":  0.0, "cdual": 0.0, "cprimal": 0.0, "obj": 0.0, "obj_": 0.0}


def add(stat, entry):
    stat["solved"] += 1 if  float(entry["absgap"]) < 1.000001 else 0
    stat["solution"] += not entry["missing"]
    stat["total"] += 1
    if not (entry["gap"] >= 0 and entry["gap"] <= 100):
        print(entry["gap"])
    assert( entry["gap"] >=0 and entry["gap"] <= 100)
    for key in sgm_keys:
        if key in entry:
            stat[key] += np.log(float(entry[key])+ shift[key])

def parsed(name):
    return name.replace("_", "\_")


def avgStat(stat):
    for key in sgm_keys:
        if stat["total"] > 0:
            stat[key] = np.exp(stat[key] / stat["total"])- shift[key]
        else:
            stat[key] = defualt_entry[key]
        val = stat[key]
        if key in int_keys:
            val = int(val)
        else:
            val = round(val,1)
        stat[key] = val


    display = ""

    for key in display_keys:
        val = stat[key]
        if key in int_keys:
            val = int(stat[key])
        else:
            val = round(val,1)
        display += str(key) + ":" +  " & "


    #print(avg_info , solver)
    #print(display)
    return display

def gettab(stat, end):
    bench_tab = ""
    bench_tab += str(round(stat["time"],1)) + " & "
    bench_tab += str(round(stat["gap"],1)) +  str("\%")  + " & "
    bench_tab += str(round(stat["obj"],1)) + str("\%")  + " & "
    bench_tab +=  "\\multicolumn{1}{l|}{" + str(stat["solved"])  + "/"  + str(stat["solution"])  +"}"
    bench_tab += " \\\\ \n" if end else " & "
    return bench_tab


def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (s[n//2-1]/2.0+s[n//2]/2.0, s[n//2])[n % 2] if n else None

def avg(lst):
    return sum(lst)/len(lst)

def printtable(algorithms_, benchmarks, has_obj_):

    detail = ""

    benchdict = {}

    bench_instances = {}
    bench_entries = {}

    instances_entries = {}

    instances_stats = {}

    path = Path(os.getcwd())

    for benchmark in benchmarks:
        # parse all results and logs, create solution for each instance
        bench_dir_path = str(path.parent.absolute()) + "/benchmarks/" + benchmark
        instances = os.listdir(bench_dir_path)
        bench_instances[benchmark] = instances

        result_dir_path = str(path.parent.absolute()) + "/results/" + benchmark
        results = os.listdir(result_dir_path)

        bench_entries[benchmark] = []
        work_entries = []
        for result in results:
            if result.find(".None.") != -1:
                entry = extractInstance(result_dir_path + "/" + result)
                entry["coverage"] =  "Small" if result.find("Small") != -1 else "Large"
                instances_entries[(entry["instance"],  entry["coverage"])] = entry
            else:
                entry = extractInstanceResult(result_dir_path + "/" + result)
                entry["coverage"] =  "Small" if result.find("Small") != - 1 else "Large"
                entry["type"] = 1
                work_entries.append(entry)

        # fill non solution and normalize
        for instance in instances:
            for cover in coverages:
                if cover == "Small":
                    detail +=  "\\multirow{2}{*}{\\texttt{"+ str(instance).replace(".txt", "").replace("_", "") +"}}"
                detail +=  "&" +  ("S" if cover == "Small" else "L") + " & " + str(int(instances_entries[(instance, cover)]["sdb_node"])) + " & " + str(instances_entries[(instance, cover)]["sdb_edge"])
                for algo in algorithms_:
                    is_find = False
                    for entry in work_entries:
                        if entry["instance"] == instance and entry["coverage"] == cover and entry["formulation"] == algo:
                            is_find = True
                            val = entry["obj"]
                            if entry["missing"]:
                                entry["obj"] = 100
                            else:
                                entry["obj"] = val / instances_entries[(instance, cover)]["sdb_node"] * 100
                            entry["isnotfind"] = False
                            bench_entries[benchmark].append(entry)
                            detail += "&" + str(round(entry["time"]/1800,1)) + " & "  + str(round(entry["gap"],1)) + "\% & " +  str(round(entry["obj"],1)) + "\%"
                    if not is_find:
                        entry = copy.copy(defualt_entry)
                        entry["formulation"]  = algo
                        entry["instance"] = instance
                        entry["coverage"] = cover
                        entry["obj"] = 100
                        entry["gap"] = 100
                        entry["missing"] = True
                        entry["isnotfind"] = True
                        bench_entries[benchmark].append(entry)
                        detail += "&" + str("-") + " & "  + str("-") + " & " +  str("-")
                detail += "\\\\" + "\n"

    file = open("detail.txt", "w")
    file.write(detail)
    file.close()



    bench_instances = {}
    bench_instances["Small"] = []
    bench_instances["Large"] = []
    minlarge = 10000000
    maxsmall = -1
    edge_list = {}
    edge_list["Small"] = []
    edge_list["Large"] = []
    dense_list = {}
    dense_list["Small"] = []
    dense_list["Large"] = []
    for k in instances_entries:
        entry = instances_entries[k]
        if entry["org_edge"] <= 150:
            maxsmall= max(maxsmall, entry["org_edge"])
            bench_instances["Small"].append(entry["instance"])
            edge_list["Small"].append(entry["org_edge"])
            dense_list["Small"].append( 2.0 * entry["org_edge"] / entry["org_node"] * (entry["org_node"] - 1))
        else:
            minlarge= min(minlarge, entry["org_edge"])
            bench_instances["Large"].append(entry["instance"])
            edge_list["Large"].append(entry["org_edge"])
            dense_list["Large"].append( 2.0 * entry["org_edge"] / entry["org_node"] * (entry["org_node"] - 1))

    print(len(bench_instances["Small"]), " ", len(bench_instances["Large"]))
    print("Small bench: min/med/max/avg/avg_dense:", min(edge_list["Small"] ), " ", median(edge_list["Small"] ), "  ", max(edge_list["Small"] ), "  ", avg(edge_list["Small"] ),"  ", avg(dense_list["Small"] ))
    print("Large bench: min/med/max/avg/avg_dense:", min(edge_list["Large"] ), " ", median(edge_list["Large"] ), "  ", max(edge_list["Large"] ), "  ", avg(edge_list["Large"] ),"  ", avg(dense_list["Large"] ))

    bench_instances["Small"] = set(bench_instances["Small"])
    bench_instances["Large"] = set(bench_instances["Large"])


    bench_entries["Small"] = []
    bench_entries["Large"] = []
    for benchmark in benchmarks:
        for entry in bench_entries[benchmark]:
            if entry["instance"] in bench_instances["Small"]:
                bench_entries["Small"].append(entry)
            else:
                bench_entries["Large"].append(entry)

    test_benchs = ["Small", "Large"]
    for benchmark in test_benchs:
        # benchmark_wise statistics
        for cover in coverages:
            for algo in algorithms_:
                bench_tab = ""
                stat = Stat(algo, cover)
                #print(algo, cover)
                for entry in bench_entries[benchmark]:
                    if entry["formulation"] == algo and entry["coverage"] == cover:
                        #print(entry)
                        add(stat, entry)
                        #print(entry)
                #print(stat)
                display = avgStat(stat)
                #print(display)
                benchdict[(benchmark, cover, algo)] = gettab(stat, False if cover == "Small" else True)
                #bench_tab = gettab(bench_tab, stat, has_obj_ , 0)
            #print("\n")
            #bench_tab += "\n"

    for benchmark, subbench_name in zip(test_benchs, ["Small", "Large"]):
        tab = ""
        for algo_ in algorithms_:
            tab += algonm_map[algo_] + " & "
            for cover_ in coverages:
                tab += benchdict[(subbench_name, cover_, algo_)]
        print(tab, "\n")

    for benchmark, subbench_name in zip(test_benchs, ["Small", "Large"]):
        for cover_ in coverages:
            ys = {}
            for algo_ in comp_algos:
                ys[algo_] = []
            for algo_ in comp_algos:
                for entry in bench_entries[benchmark]:
                    if entry["formulation"] == algo_ and entry["coverage"] == cover_:
                        ys[algo_].append( entry["gap"] )
            for algo_, style, color, marker in zip(comp_algos, styles, colors, markers):
                #print(ys[algo_])
                ys[algo_].sort()
                j = 0
                for y in  ys[algo_]:
                    if y < 1e-1:
                        j += 1
                k = len(ys[algo_])
                xs = [i + j for i in range(k - j)]
                print(len(ys[algo_][j:]), len(xs))
                plt.plot(ys[algo_][j:], xs, label = algonm_map[algo_], linestyle =style, color = color, marker = marker,  markersize=4)
            #plt.ylim(ymin=100)
            plt.legend(loc='lower right')
            plt.axis('tight')
            plt.xlabel('relative gap %')
            plt.ylabel('#instances')
            plt.savefig(benchmark + cover_ + ".pdf", format="pdf", bbox_inches="tight")
            plt.show()


# display table 2
print("table \n")
algorithms_  = [ algo for algo in algorithms if algo != "None"]
has_obj_ = False
printtable(algorithms_, benchmarks, has_obj_)

