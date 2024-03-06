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
time_limit =1800
algorithms  = ["EF", "EVF", "EFP0", "EFP", "EFPV", "EFPV2", "EFPLg", "EFPD", "EVFP0", "EVFP", "EVFPV", "EVFPL", "None"]
coverages = ["Small", "Large"]
benchmarks = ["city", "Kgroup_A", "Kgroup_B", "random_A", "random_B"]


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
    stat_keys = ["time", "node", "gap", "obj", "bound", "instance", "algo"]
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
    entries["algo"] = stat_dict["algo"].split()[1]
    entries["missing"] = False
    return entries
    #writer.writerow([name, algo, primal_bound, gap, total_time, pricing_time, pricing_calls, columns_gen, nodes])

def extractInstance(file_path):
    ls = open(file_path).readlines()
    #print(file)
    stat_keys = ["instance", "algo", "dlt", "org_node", "org_edge", "org_min_len", "org_avg_len", "org_max_len",  "sdb_node", "sdb_edge", "sdb_min_len", "sdb_avg_len", "sdb_max_len",  "dtf_node", "dtf_edge", "dtf_min_len", "dtf_avg_len", "dtf_max_len"]
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
    entries["algo"] = stat_dict["algo"].split()[1]
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
    stat["solution"] += entry["type"]
    stat["total"] += 1
    if not (entry["gap"] >= 0  and entry["gap"] <= 100):
        pass#print(entry["gap"])
    assert( entry["gap"] >=0  and entry["gap"] <= 100)
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

def gettab( stat,  has_obj_, end):
    bench_tab = ""
    bench_tab += str(round(stat["time"],1)) + " & "
    bench_tab += str(round(stat["gap"],1)) +  str("\%")  + " & "
    bench_tab += str(round(stat["obj"],1)) + str("\%")  + " & "
    if has_obj_:
        bench_tab += str(round(stat["obj_"],1)) + str("\%") + " & "
    if not end:
        bench_tab +=   str(stat["solved"])  + "/"  + str(stat["solution"])  + "/" + str(stat["total"])  + " & "
    else:
        bench_tab +=   "\\multicolumn{1}{l|}{" + str(stat["solved"])  + "/"  + str(stat["solution"])  + "/" + str(stat["total"]) +"}"
    return bench_tab


def addaxes(axes, i, algo1, algo2, dual_results, primal_results):
    #print(dual_results[algo1], dual_results[algo2])
    axes[i,0].scatter(dual_results[algo1],dual_results[algo2], color = 'blue', marker = '+')
    axes[i,0].plot(dual_results[algo1],dual_results[algo1], color = 'green')
    axes[i,0].set_xlabel(algo1)
    axes[i,0].set_ylabel(algo2)
    axes[i,0].set_title("relative dual gap")

    axes[i,1].scatter(primal_results[algo1],primal_results[algo2], color = 'red', marker = 'x')
    axes[i,1].plot(primal_results[algo1],primal_results[algo1], color = 'orange')
    axes[i,1].set_xlabel(algo1)
    axes[i,1].set_ylabel(algo2)
    axes[i,1].set_title("relative primal bound")

def addaxes_(axes, i, algo1, algo2, dual_results, primal_results):
    axes[0].scatter(dual_results[algo1],dual_results[algo2], color = 'blue', marker = '+')
    axes[0].plot(dual_results[algo1],dual_results[algo1], color = 'green')
    axes[0].set_xlabel(algo1)
    axes[0].set_ylabel(algo2)
    axes[0].set_title("relative dual gap")

    axes[1].scatter(primal_results[algo1],primal_results[algo2], color = 'red', marker = 'x')
    axes[1].plot(primal_results[algo1],primal_results[algo1], color = 'orange')
    axes[1].set_xlabel(algo1)
    axes[1].set_ylabel(algo2)
    axes[1].set_title("relative primal bound")

def printtable(algorithms_, has_obj_):
    normal_entries = []
    extreme_entries = []

    details = ""

    benchdict = {}  

    for benchmark in benchmarks:
        # parse all results and logs, create solution for each instance
        bench_tab = ""
        path = Path(os.getcwd())
        bench_dir_path = str(path.parent.absolute()) + "/results/" + benchmark
        #print(bench_dir_path)
        instances = os.listdir(bench_dir_path)
        entries = []
        instance_stats = []
        for instance in instances:
            #print(instance, "Small" if instance.find("Small") != -1 else "Large")
            if instance.find(".None.") != -1:
                entry = extractInstance(bench_dir_path + "/"+instance)
                entry["coverage"] =  "Small" if instance.find("Small") != -1 else "Large"
                entry["dual"] = max_dual_bound
                entry["primal"] = min_primal_bound
                instance_stats.append(entry)
            else:
                entry = extractInstanceResult(bench_dir_path + "/"+instance)
                entry["coverage"] =  "Small" if instance.find("Small") != - 1 else "Large"
                entry["type"] = 1
                entries.append(entry)


        # compute over-all stat
        for instance_stat in instance_stats:
            for entry in entries:
                if entry["instance"] == instance_stat["instance"] and entry["coverage"] == instance_stat["coverage"]:
                    instance_stat["dual"] = max(entry["bound"], instance_stat["dual"])
                    instance_stat["primal"] = min(entry["obj"], instance_stat["primal"])


        # compute closed dual and primal
        for instance_stat in instance_stats:
            for algo in algorithms_:
                for entry in entries:
                    if entry["instance"] == instance_stat["instance"] and entry["coverage"] == instance_stat["coverage"]:
                        entry["cdual"] = max( abs(instance_stat["dual"]  - instance_stat["primal"]) ,1e-6) / max( abs(entry["bound"]  - instance_stat["primal"]) ,1e-6)
                        entry["cprimal"] = max( abs(instance_stat["dual"]  - instance_stat["primal"]) ,1e-6) / max( abs(instance_stat["dual"]  - entry["obj"]) ,1e-6)
            

        #print(entries)
        # fill non solution
        for instance_stat in instance_stats:
            for algo in algorithms_:
                is_find = False
                for entry in entries:
                    if entry["instance"] == instance_stat["instance"] and entry["coverage"] == instance_stat["coverage"] and entry["algo"] == algo:
                        is_find = True
                        val = entry["obj"] 
                        entry["obj"] = val / instance_stat["org_node"] * 100
                        entry["isnotfind"] = False
                if not is_find:
                    entry = copy.copy(defualt_entry)
                    entry["algo"]  = algo
                    entry["instance"] = instance_stat["instance"]
                    entry["coverage"] = instance_stat["coverage"]
                    entry["obj"] = 100
                    entry["isnotfind"] = True
                    entries.append(entry)
                
        if benchmark == "Kgroup_B":
            extreme_entries += entries
        else:
            normal_entries += entries



        # benchmark_wise statistics 
        for cover in coverages:
            for algo in algorithms_:
                bench_tab = ""
                stat = Stat(algo, cover)
                #print(algo, cover)
                for entry in entries:
                    if entry["algo"] == algo and entry["coverage"] == cover:
                        add(stat, entry)
                        #print(entry)
                #print(stat)
                display = avgStat(stat)
                #print(display)
                #bench_tab = gettab(bench_tab, stat, has_obj_ , 0)
                benchdict[(benchmark, cover, algo ,0 )] = gettab( stat, has_obj_ , 0)
                benchdict[(benchmark, cover, algo ,1 )] = gettab( stat, has_obj_ , 1)
            #print("\n")
            #bench_tab += "\n"

        #bench_tab = bench_tab[0:-2]
        #print(bench_tab, "\n\n")

                
        #for cover in coverages:
        #    max_node = 0
        #    min_node = 2**30
        t = 0
        for instance_stat_ in instance_stats:
            if instance_stat_["coverage"] == "Small":
                continue
            name = str(parse_name(instance_stat_["instance"])[0:-4])
            detail =  "\\multirow{2}{*}{\\texttt{"+ str(parse_name(instance_stat_["instance"])[0:-4]) +"}}"
            for cover in coverages:
                has_cover = False
                for instance_stat in instance_stats:
                    if instance_stat["coverage"] == cover and instance_stat_["instance"][0:-4] == instance_stat["instance"][0:-4]:
                        detail +=  "&" + cover + " & " + str(int(instance_stat["sdb_node"])) + " & " + str(int(instance_stat["sdb_edge"])) + " & " + str(int(instance_stat["dtf_node"])) + " & " + str(int(instance_stat["dtf_edge"]))
                        for algo in algorithms_:
                                if algo == "None":
                                    continue
                                iffind = False
                                for entry in entries:
                                    if entry["instance"] ==  instance_stat["instance"] and entry["algo"] == algo and entry["coverage"] == cover and "isnotfind" in entry and not entry["isnotfind"]:
                                        detail += "&" + str(round(entry["time"],1)) + " & "  + str(round(entry["gap"],1)) + "\% & " +  str(round(entry["obj"],1)) + "\%"
                                        iffind = True
                                if not iffind:
                                        detail += "&" + str("-") + " & "  + str("-") + " & " +  str("-") 
                        detail += "\\\\" + "\n"
                        has_cover = True
                #print(cover, has_cover, detail)
            #t +=  1
            #print(detail)
            details += detail
                        #max_node = max(instance_stat["org_node"], max_node)
                        #min_node = min(instance_stat["org_node"], min_node)
                
            #print(max_node, min_node)

    alldisentries = []
    for cover_ in coverages:
        for algo_ in algorithms_:
            stat = Stat(algo_, cover_)
            for benchmark in benchmarks:
                # parse all results and logs, create solution for each instance
                path = Path(os.getcwd())
                bench_dir_path = str(path.parent.absolute()) + "/results/" + benchmark
                #print(bench_dir_path)
                instances = os.listdir(bench_dir_path)
                entries = []
                instance_stats = []
                for instance in instances:
                    #print(instance, instance.find(".None."))
                    if instance.find(".None.") != -1:
                        entry = extractInstance(bench_dir_path + "/"+instance)
                        entry["coverage"] =  "Small" if instance.find("Small") != -1 else "Large"
                        entry["dual"] = max_dual_bound
                        entry["primal"] = min_primal_bound
                        instance_stats.append(entry)
                    else:
                        entry = extractInstanceResult(bench_dir_path + "/"+instance)
                        entry["coverage"] =  "Small" if instance.find("Small") != -1 else "Large"
                        entry["type"] = 1
                        entries.append(entry)

                # compute over-all stat
                for instance_stat in instance_stats:
                    for entry in entries:
                        if entry["instance"] == instance_stat["instance"] and entry["coverage"] == instance_stat["coverage"]:
                            instance_stat["dual"] = max(entry["bound"], instance_stat["dual"])
                            instance_stat["primal"] = min(entry["obj"], instance_stat["primal"])
                
                # compute over-all stat
                for instance_stat in instance_stats:
                    for entry in entries:
                        if entry["instance"] == instance_stat["instance"] and entry["coverage"] == instance_stat["coverage"]:
                            instance_stat["dual"] = max(entry["bound"], instance_stat["dual"])
                            instance_stat["primal"] = min(entry["obj"], instance_stat["primal"])

                # compute closed dual and primal
                for instance_stat in instance_stats:
                    for algo in algorithms_:
                        for entry in entries:
                            if entry["instance"] == instance_stat["instance"] and entry["coverage"] == instance_stat["coverage"]:
                                entry["cdual"] = max( abs(instance_stat["dual"]  - instance_stat["primal"]) ,1e-6) / max( abs(entry["bound"]  - instance_stat["primal"]) ,1e-6)
                                entry["cprimal"] = max( abs(instance_stat["dual"]  - instance_stat["primal"]) ,1e-6) / max( abs(instance_stat["dual"]  - entry["obj"]) ,1e-6)
                                
                # fill non solution
                change_mode = False
                if "SFD" in algorithms_ and  "RF" in algorithms_ and len(algorithms_) == 2:
                    change_mode = True
                    for instance_stat in instance_stats:
                        for entry in entries:
                            if entry["instance"] == instance_stat["instance"] and entry["coverage"] == instance_stat["coverage"] and entry["algo"] == "SFD":
                                instance_stat["norm"] = entry["obj"] 
                for instance_stat in instance_stats:
                    for algo in algorithms_:
                        is_find = False
                        for entry in entries:
                            if entry["instance"] == instance_stat["instance"] and entry["coverage"] == instance_stat["coverage"] and entry["algo"] == algo:
                                is_find = True
                                if change_mode:
                                    val = entry["obj"] 
                                    entry["obj"] = val / instance_stat["org_node"] * 100
                                    entry["obj_"] = val / instance_stat["norm"] * 100
                                else:
                                    val = entry["obj"] 
                                    entry["obj"] = val / instance_stat["org_node"] * 100
                            entry["isnotfind"] = False
                        if not is_find:
                            entry = copy.copy(defualt_entry)
                            entry["algo"]  = algo
                            entry["instance"] = instance_stat["instance"]
                            entry["coverage"] = instance_stat["coverage"]
                            entry["obj"] = 100
                            entry["isnotfind"] = True
                            entries.append(entry)
                alldisentries += entries

                for entry in entries:
                    if entry["algo"] == algo_ and entry["coverage"] == cover_:
                        add(stat, entry)
            display = avgStat(stat)
            #print(stat)
            benchdict[("all", cover_, algo_, 0 )] = gettab(stat, has_obj_ , 0)
            benchdict[("all", cover_, algo_, 1 )] = gettab(stat, has_obj_ , 1)
            #allbench_tab += "\n"
    #print(allbench_tab)

    algos = copy.copy(algorithms_)
    if len(algos) % 2 == 0:
        nrows_ = len(algos) // 2
    else:
        nrows_ = len(algos) // 2
    fig, axes = plt.subplots(nrows=nrows_, ncols=2, figsize=(10, 6))  # define the figure and subplots
    pairs = []

    print(algos, algorithms_)

    for i in range(nrows_):
        pairs.append((algos[2*i], algos[2*i + 1]))

    i = 0
    for pair in pairs:
        algo1 = pair[0]
        algo2 = pair[1]
        dual_results = {algo1:[], algo2:[]}
        primal_results = {algo1:[], algo2:[]}
        instances = []
        for entry in normal_entries:
            if entry["algo"] == algo1 and not entry["isnotfind"]:
                has_entry = False
                for entry_ in normal_entries:
                    if entry_["algo"] == algo2 and entry_["instance"] == entry["instance"] and not entry["isnotfind"]:
                        has_entry = True
                        dual_results[algo1].append(entry["gap"] / 100.0)
                        dual_results[algo2].append(entry_["gap"] / 100.0)
                        primal_results[algo1].append(entry["obj"] / 100.0)
                        primal_results[algo2].append(entry_["obj"] / 100.0)
        addaxes(axes, i, algo1, algo2, dual_results, primal_results)
        i+=1
    fig.tight_layout()
    plt.savefig('scatter5.pdf') 

    benchmarks_ = benchmarks + ["all"]
    groudsize = int(np.floor(len(algorithms_) / 2))
    for groudid in range(groudsize):
        algos = [algorithms_[2*groudid], algorithms_[2*groudid + 1]]
        print(groudid, " ", algos,"\n")
        tab = ""
        for benchmark in benchmarks_:
            tab+= "\\multirow{2}{*}{\\texttt{"+ parsed(benchmark) +"}}"
            for cover in  coverages:
                tab += "&" + cover + " & " 
                for algo in algos:
                    if algo == algos[-1]:
                    #print(benchdict[(benchmark, cover, algo)], "\n")
                        s = benchdict[(benchmark, cover, algo, 1)]
                        tab += s
                    else:
                        tab += benchdict[(benchmark, cover, algo, 0)]
                tab += " \\\\"
                tab += "\n"
        print(tab)
        print("\n")  
    if len(algorithms_)- 2 *(groudsize) != 0:
        tab = ""
        algos = [algorithms_[-1]]
        print( " ", algos,"\n")
        for benchmark in benchmarks_:
            tab+= "\\multirow{2}{*}{\\texttt{"+ parsed(benchmark) +"}}"
            for cover in  coverages:
                tab += "&" + cover + " & " 
                for algo in algos:
                    if algo == algos[-1]:
                    #print(benchdict[(benchmark, cover, algo)], "\n")
                        s = benchdict[(benchmark, cover, algo, 1)]
                        tab += s
                    else:
                        tab += benchdict[(benchmark, cover, algo, 0)]
                tab += " \\\\"
                tab += "\n"
        print(tab)
    return details

# display table 2
print("table 2\n")
algorithms_  = [ algo for algo in algorithms if algo != "None"]
has_obj_ = False
printtable(algorithms_, has_obj_)


# write details

algorithms_  = [ algo for algo in algorithms if algo != "None"]
has_obj_ = False
details = printtable(algorithms_, has_obj_)

file = open("details.txt", "w")
file.write(details)
file.close()
