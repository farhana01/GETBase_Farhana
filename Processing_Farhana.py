import math
import os
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np
# Need to install on linux before use, and include as packages in PyCharm
from graphviz import Graph  # Network visualizer (nodes and edges)

import copy

raw_data_path = "./Output/"
figure_path = "./Figures/"
main_filename = "best.dat"
num_runs = 30
lower_better = False
precision = 6
col_width = 6 + precision
num_nodes = 256


# Go through a best##.dat file and gather all important results from the file.
def get_data(file_path: str):
    data = []
    fits = []
    hists = []
    nets_edge_lists = []
    epi_prof_flag = False
    graph_flag = False
    edge_lists_flag = False
    variant_parents = []
    nets_num_edges = []
    # variant_parents[run_num][variant_idx] = idx of the variant that gave birth to the variant with index variant_idx in run run_num

    variant_starts = []
    variant_profs = []

    with open(file_path) as f:
        lines = f.readlines()
        for line in lines:
            if line.__contains__("Best Fitness: "):
                line = line.split(": ")[1]  # i.e. line = "20387.8\n"
                # Removes all the whitespace and newline/carriage return from the right side of the string.
                line = line.rstrip()  # i.e. line = "20387.8"
                fits.append(float(line))
                pass
            elif line.__contains__("0.5"):
                epi_prof_flag = False
                pass
            elif line.__contains__("Graph"):
                graph_flag = True
                pass
            elif line.__contains__("W Hist:"):
                graph_flag = False
                edge_lists_flag = True
                one_net_lists = []
                pass
            elif line.rstrip() == "":
                edge_lists_flag = False
                nets_edge_lists.append(one_net_lists)
                pass
            elif line.__contains__("Severity Histogram: "):
                line = line.split(": ")[1]
                line = line.rstrip()
                line = line.split("\t")
                one_hist = [int(x) for x in line]
                hists.append(one_hist)
                epi_prof_flag = False
                variant_parents.append(one_run_var_par)
                variant_starts.append(one_run_var_starts)
                variant_profs.append(one_run_var_profs)
                pass
            elif line.__contains__("Epidemic Profile"):
                epi_prof_flag = True
                one_run_var_par = []
                one_run_var_starts = []
                one_run_var_profs = []
                pass
            elif epi_prof_flag:
                if line.__contains__("->"):
                    line = line.split("->")
                    var_par = line[0].replace("-", "").rstrip()
                    line = line[1].split("[")[1].split("-")
                    start_day = line[0].rstrip()
                    line = line[1].split("]:")[1].strip().split("\t")
                    var_prof = [int(x) for x in line]
                    # Get Variant Parent's
                    if var_par == "NA":
                        one_run_var_par.append(-1)
                        pass
                    else:
                        var_par = var_par.split("V")[1]
                        one_run_var_par.append(int(var_par))
                        pass

                    one_run_var_starts.append(int(start_day))
                    one_run_var_profs.append(var_prof)
                    pass
                pass
            elif graph_flag:
                if line.__contains__("Edges: "):
                    nets_num_edges.append(int(line.split(": ")[1].rstrip()))
                    pass
                pass
            elif edge_lists_flag:
                one_net_lists.append([int(i) for i in line.rstrip().split(" ")])
                pass
            pass
        pass

    # data[run_num][0] = run number
    # data[run_num][1] = run's fitness
    # data[run_num][2][:] = run's best epidemic infection severity histogram
    # data[run_num][3][:] = variants' parent's index for all variants in the best epidemic for this run
    # data[run_num][4][:] = variants' start date for all variants in the best epidemic for this run
    # data[run_num][5][:][:] = epidemic profiles for all variants in the best epidemic for this run
    # data[run_num][6] = number of edges in the best network for this run
    # data[run_num][7][:][:] = edge lists for the best network for this run
    data = [[idx + 1, fits[idx], hists[idx], variant_parents[idx], variant_starts[idx], variant_profs[idx],
             nets_num_edges[idx], nets_edge_lists[idx]] for idx in range(num_runs)]
    data.sort(key=itemgetter(1))  # Ascending
    data.reverse()  # Descending

    return data


def make_histogram(vals: [], exp_name: str):
    f = plt.figure()
    f.set_dpi(600)
    f.set_figheight(3)
    f.set_figwidth(8)
    plot = f.add_subplot(111)

    x_locs = [idx + 1 for idx in range(len(vals[1:35]))]
    x_lbls = [str(i) for i in x_locs]
    p = plot.bar(x_locs, vals[1:35], label=x_lbls)
    plot.bar_label(p, label_type='edge', fontsize=8)

    y_max = plot.get_ylim()[1]
    x_max = plot.get_xlim()[1]
    plot.set_ylim(0, y_max)
    plot.set_xlim(0, x_max)

    plot.set_title("Best Epidemic Histogram Using " + exp_name, fontsize=12)
    plot.set_xlabel("Severity of Infection", fontsize=8)
    plot.set_ylabel('Number of Infections', fontsize=8)
    plot.minorticks_on()
    plot.grid(visible="True", axis="x", which='minor', color="lightgray", linewidth=0.5)
    plot.grid(visible="True", axis="x", which='major', color="white", linewidth=0.5)
    plot.grid(visible="True", axis="y", which='major', color="white", linewidth=0.5)
    f.subplots_adjust(bottom=0.12, top=0.91, left=0.07, right=0.99)
    f.savefig(figure_path + exp_name + "_hist.png", dpi=600)
    plt.close()
    pass


def make_epi_profile(variant_parents: [], variant_starts: [], variant_profs: [], exp_name: str):
    f = plt.figure()
    f.set_dpi(600)
    f.set_figheight(4)
    f.set_figwidth(8)

    plot = f.add_subplot(111)

    var_xs = []
    var_ys = []

    for idx in range(len(variant_parents)):
        var_xs.append([i for i in range(variant_starts[idx], variant_starts[idx] + len(variant_profs[idx]) + 1)])
        var_prof = variant_profs[idx]
        var_prof.append(0)
        var_ys.append(var_prof)
        pass

    # Legend
    for idx in range(len(variant_parents)):
        if idx == 0:
            lbl = "V" + str(idx) + " (NA)"
            pass
        else:
            lbl = "V" + str(idx) + " (V" + str(variant_parents[idx]) + ")"
            pass
        plot.plot(var_xs[idx], var_ys[idx], label=lbl)
        pass

    y_max = plot.get_ylim()[1]
    x_max = plot.get_xlim()[1]
    plot.set_ylim(0, y_max)
    plot.set_xlim(0, x_max)

    plot.set_title("Best Epidemic Curve Using " + exp_name, fontsize=12)
    plot.set_xlabel("Day", fontsize=8)
    plot.set_ylabel("New Infections", fontsize=8)
    plot.minorticks_on()
    plot.grid(visible="True", axis="x", which='minor', color="darkgray", linewidth=0.5)
    plot.grid(visible="True", axis="x", which='major', color="black", linewidth=0.5)
    plot.grid(visible="True", axis="y", which='major', color="black", linewidth=0.5)
    plt.legend(title="Variant ID (Parent)", borderaxespad=0, bbox_to_anchor=(1, 0.5),
               loc="center left", fontsize=6, title_fontsize=6, ncol=2)
    f.subplots_adjust(bottom=0.085, top=0.93, left=0.06, right=0.8)
    f.savefig(figure_path + exp_name + "_profile.png", dpi=600)
    plt.close()
    pass


def make_lines(the_ys: [], labels: [], title: str):
    f = plt.figure()
    f.set_dpi(600)
    f.set_figheight(4)
    f.set_figwidth(8)

    plot = f.add_subplot(111)

    xs = [i for i in range(1, 36)]

    for idx in range(len(the_ys)):
        plot.plot(xs, the_ys[idx][1:36], label=labels[idx])
        pass

    y_max = plot.get_ylim()[1]
    x_max = plot.get_xlim()[1]
    plot.set_ylim(0, y_max)
    plot.set_xlim(0, x_max)

    plot.set_title(title, fontsize=12)
    plot.set_xlabel("Severity of Infection", fontsize=8)
    plot.set_ylabel('Number of Infections', fontsize=8)
    plot.minorticks_on()
    plot.grid(visible="True", axis="x", which='minor', color="lightgray", linewidth=0.5)
    plot.grid(visible="True", axis="x", which='major', color="white", linewidth=0.5)
    plot.grid(visible="True", axis="y", which='major', color="white", linewidth=0.5)
    plt.legend(title="Legend", loc="best", fontsize=6, title_fontsize=6)
    f.subplots_adjust(bottom=0.08, top=0.93, left=0.07, right=0.99)
    f.savefig(figure_path + title + ".png", dpi=600)
    plt.close()
    pass


def calc_avg_inf_var(data: [], info: str):
    total_inf = 0
    total_days = 0
    total_vars = 0
    for dat in data:
        total_vars += len(dat[5])
        for var in dat[5]:
            total_inf += sum(var)
            total_days += len(var)
            pass
        pass
    print("Average Inf/Day/Var for " + info + ": " + str(float(total_inf / total_days)))
    print("Average Inf/Var for " + info + ": " + str(float(total_inf / total_vars)))
    print("Average Total Infections In Epidemic for " + info + ": " + str(float(total_inf / num_runs)))
    print("Average Variant Length for " + info + ": " + str(float(total_days / total_vars)))
    pass


def edge_list(node_lists, num_nodes):
    edge_lists = []
    for fr in range(num_nodes):
        for to in node_lists[fr]:
            edge_lists.append([fr, to])
            pass
        pass
    return edge_lists


def make_graph(el: [], out_file: str, verts: int):
    g = Graph(engine="sfdp")
    e_cout = 0

    g.graph_attr.update(dpi='300', size="6,6", splines='true', overlap='false')
    # g.node_attr.update(color='black', shape='square', width='0.02', height='0.02')
    g.node_attr.update(color='black', shape='point', fixedsize='true', width='0.01')
    g.edge_attr.update(color='black', penwidth='0.25')

    for n in range(verts):
        if n == 0:
            g.node(str(n), color='red')
            pass
        else:
            g.node(str(n))
        pass

    for idx, d in enumerate(el):
        if d[0] < d[1]:
            g.edge(str(d[0]), str(d[1]), color='black')
            e_cout += 1
            pass
        pass
    g.render(filename=out_file, directory=figure_path, cleanup=True, format='png')
    del g
    print("Made network: " + out_file + " with " + str(e_cout) + " Edges")
    pass


def node_deg_hist(el: []):
    degs = []
    for li in el:
        degs.append(len(li))
        pass

    plt.hist(degs)
    plt.show()

    pass


def make_box(data, labels, title, out_path, vert_lbl):
    colours = ["#02a5ff", "#ff3300", "#66cc33", "#ffcc00", "#123456"]

    data = copy.deepcopy(data)
    plt.style.use("seaborn-v0_8")
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    # Create a figure object
    f = plt.figure()
    f.set_figheight(8)
    f.set_figwidth(6)
    plot = f.add_subplot(111)

    # Calculate x values and generate the plot.
    xs = [i + 1 for i in range(len(data))]
    just_fits = [[data[exp][x][1] for x in range(len(data[exp]))] for exp in range(len(data))]
    vp = plot.violinplot(just_fits, xs, showmedians=True, widths=0.85)

    # Make it fancy!
    for pc in vp["bodies"]:
        pc.set_facecolor("#5770DB")
        pc.set_linewidth(1)
        pc.set_edgecolor("black")
        pass

    # for pc, color in zip(vp['bodies'], colours):
    #     pc.set_facecolor(color)
    #     pc.set_linewidth(1)
    #     pc.set_alpha(1)
    #     pc.set_color(color)
    # for idx, col in enumerate(colours):
    #     for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
    #         vc = vp[partname][idx]
    #         vc.set_linewidth(1)
    #         vc.set_alpha(1)
    #         vc.set_color(col)
    #         pass
    #     pass

    plot.grid(visible="True", axis="y", which='minor', color="white", linewidth=0.5)
    plt.minorticks_on()
    # plt.yscale('log')

    f.suptitle(title, fontsize=14)
    plot.set_xlabel(vert_lbl, fontsize=12)
    plot.set_ylabel("Fitness", fontsize=12)

    # plt.yticks(ticks)
    plt.xticks([x + 1 for x in range(len(data))], labels=labels, rotation=25, fontsize=10)

    # f.subplots_adjust(bottom=0.25, top=0.93, left=0.1, right=0.99)
    f.tight_layout()
    f.savefig(out_path + "Boxplot " + title + ".png", dpi=500)

    plt.close()
    pass


def main():
    print("START")
    folder_names = os.listdir(raw_data_path)
    if len(folder_names) == 2:
        exp_titles = ["Epidemic Severity Fitness", "Epidemic Spread Fitness"]
        pass
    else:
        exp_titles = ["Epidemic Spread Fitness"]
        pass

    data = []
    for fold in folder_names:
        dat = get_data(raw_data_path + fold + "/" + main_filename)
        data.append(dat)
        pass

    make_box(data, [str(x) for x in range(len(data))], "Compare Results", figure_path, "Experiment")

    # data[exp_num][run_num][1] = run's fitness
    # data[exp_num][run_num][2][:] = run's histogram
    # data[exp_num][run_num][3][:] = variant's parent's index for all variants

    plt.style.use("seaborn-v0_8-dark")
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)

    avg_hists = []
    avg_edges = []

    for exp_num, title in enumerate(exp_titles):
        make_histogram(data[exp_num][0][2], title)

        avg_hist = [0 for _ in range(128)]
        total_edges = 0
        for run in range(num_runs):
            total_edges += data[exp_num][run][6]
            for idx, val in enumerate(data[exp_num][run][2]):
                avg_hist[idx] += val
                pass
            pass
        for idx in range(len(avg_hist)):
            avg_hist[idx] = int(avg_hist[idx] / num_runs)
            pass

        avg_hists.append(avg_hist)
        avg_edges.append(total_edges / num_runs)
        make_histogram(avg_hist, "AVERAGE " + title)

        make_epi_profile(data[exp_num][0][3], data[exp_num][0][4], data[exp_num][0][5], title)
        pass

    if len(folder_names) == 2:
        make_lines([avg_hists[0], avg_hists[1]], ["Severity", "Spread"], "AVERAGE Histogram for Both Fitness Functions")
        calc_avg_inf_var(data[0], "Epidemic Severity")
        calc_avg_inf_var(data[1], "Epidemic Spread")
        pass
    else:
        make_lines([avg_hists[0]], ["Spread"], "AVERAGE Histogram for Both Fitness Functions")
        calc_avg_inf_var(data[0], "Epidemic Spread")
        pass

    print("")

    if len(folder_names) == 2:
        print("Average Number of Edges for Epidemic Severity: " + str(avg_edges[0]))
        print("Average Number of Edges for Epidemic Spread: " + str(avg_edges[1]))
        node_deg_hist(data[0][0][7])
        node_deg_hist(data[1][0][7])
        make_graph(edge_list(data[0][0][7], num_nodes), "Best Network with Epidemic Severity Fitness", num_nodes)
        make_graph(edge_list(data[1][0][7], num_nodes), "Best Network with Epidemic Spread Fitness", num_nodes)
        pass
    else:
        print("Average Number of Edges for Epidemic Spread: " + str(avg_edges[0]))
        node_deg_hist(data[0][0][7])
        make_graph(edge_list(data[0][0][7], num_nodes), "Best Network with Epidemic Spread Fitness", num_nodes)
        pass

    print("END")
    pass


main()
