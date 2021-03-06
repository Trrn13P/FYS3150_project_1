import matplotlib.pyplot as plt
import numpy as np

def read_file(filename):
    infile = open(filename,"r")
    first_line = infile.readline()
    if first_line[0] == "|":
        N = eval(first_line.split()[1][2:])
        epsilon_max = eval(first_line.split()[3][12:])
        epsilon_tot = eval(first_line.split()[5][12:])
        cpu_time = eval(first_line.split()[7][9:])
        log_10_h = eval(first_line.split()[9][10:])

        x = []; v = []; u = []
        infile.readline()
        for line in infile:
            inline = line.split()
            x.append(eval(inline[0]))
            v.append(eval(inline[1]))
            u.append(eval(inline[2]))
        if filename[7:10] == "spe":
            FLOPS = 7*N
        if filename[7:10] == "gen":
            FLOPS = N*11
        infile.close()
        return x,v,u,FLOPS,N,epsilon_max,epsilon_tot,cpu_time,log_10_h

    if first_line[0] == "N":
        N = []; epsilon_max = []; log_10_h = []; cpu_time = [];
        N.append(eval(first_line.split()[0][2:]))
        epsilon_max.append(eval(first_line.split()[1][12:]))
        log_10_h.append(eval(first_line.split()[2][10:]))
        cpu_time.append(eval(first_line.split()[3][9:]))
        for line in infile:
            N.append(eval(line.split()[0][2:]))
            epsilon_max.append(eval(line.split()[1][12:]))
            log_10_h.append(eval(line.split()[2][10:]))
            cpu_time.append(eval(line.split()[3][9:]))
        return N, epsilon_max, log_10_h, cpu_time


filenames_1 = ["./data/genN10.txt","./data/genN100.txt","./data/genN1000.txt",\
"./data/speN10.txt","./data/speN100.txt","./data/speN1000.txt"]

filenames_2 = ["./data/gen_stats.txt","./data/spe_stats.txt"]


for filename in filenames_1:
    x,v,u,FLOPS,N,epsilon_max,epsilon_tot,cpu_time,log_10_h = read_file(filename)
    plt.xlabel("x")
    plt.ylabel("u(x), v(x)")

    plt.plot(x,v,label=r"Numerical solution, $n={:}$ steps".format(N))
    plt.plot(x,u,label="Analytic solution")
    if filename[7:10] == "spe":
        type = "Special algo"
    if filename[7:10] == "gen":
        type = "General algo"
    #plt.title(r"{:} $FLOPS={:}$".format(type,FLOPS))
    plt.legend() ; plt.grid()
    argv = "save"
    if argv == "plot":
        plt.show()
    if argv == "save":
        plt.savefig("./figures/1b_{:}_{:}.png".format(type[0:3],N))
    plt.clf()

n_algos = []
cpu_time_algos = []
for filename in filenames_2:
    N, epsilon_max, log_10_h, cpu_time = read_file(filename)
    cpu_time_algos.append(cpu_time)
    n_algos.append(N)
    plt.xlabel(r"$log_{10}(h)$")
    plt.ylabel(r"$\varepsilon$")
    plt.plot(log_10_h,epsilon_max)
    if filename[7:10] == "spe":
        type = "Special algo"
    if filename[7:10] == "gen":
        type = "General algo"

    plt.title("{:}".format(type))
    plt.grid(); #plt.legend()

    argv = "save"
    if argv == "plot":
        plt.show()
    if argv == "save":
        plt.savefig("./figures/1d_{:}_eps.png".format(type[0:3]))
    plt.clf()

filenames = ["./data/LU10.txt","./data/LU100.txt","./data/LU1000.txt"]
cpu_time_LU = []
for filename in filenames:
    infile = open(filename,"r")
    first_line = infile.readline().split()
    x = []; u = []; v = []
    for i in infile:
        line = i.split()
        x.append(eval(line[0]))
        v.append(eval(line[1]))
        u.append(eval(line[2]))
    cpu_time_LU.append(eval(first_line[-1][9:]))

    plt.plot(x,v,label="Numerical, n= {:}".format(filename[9:11]))
    plt.plot(x,u,label="Analytic")
    plt.xlabel("x")
    plt.ylabel("u(x), v(x)")
    plt.grid()
    plt.legend()
    plt.savefig("./figures/LU{:}.png".format(len(x)-2))
    plt.clf()

n = []
for i in range(1,8):
    n.append(int(i))
cpu_time_gen = np.log10(cpu_time_algos[0])
cpu_time_spe = np.log10(cpu_time_algos[1])
cpu_time_LU = np.log10(cpu_time_LU)

plt.plot(n,cpu_time_gen,label="General solution")
plt.plot(n,cpu_time_spe,label="Specialized solution")
plt.plot(n[0:3],cpu_time_LU,label="LU-decomposition")

plt.xlabel("n")
plt.ylabel(r"$log_10$(CPU-time) [s]")
plt.xticks(ticks = n,labels=[r"$10^1$",r"$10^2$",r"$10^3$",r"$10^4$",r"$10^5$",r"$10^6$",r"$10^7$"])
plt.legend()
plt.savefig("./figures/CPU_times")
plt.clf()
