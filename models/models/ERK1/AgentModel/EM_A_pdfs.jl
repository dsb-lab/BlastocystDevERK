using PyPlot
using PyCall
@pyimport matplotlib.patches as patch
if pwd() !== "/Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel"
    println("changing wd to /Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel")
    cd()
    cd("/Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel")
end

include("Constants/constants_mechanical.jl")
#include("Constants/constants_biochemical_RM.jl")
include("Constants/constants_biochemical_EM.jl")
include("Utils/general_utils.jl")
include("Utils/agentmodel_utils.jl")

subgroups=4
nn=10
NANOG, GATA6, FGF, tots = compute_pdfs(subgroups, nn)
NANOGmax = maximum(maximum.(NANOG))
NANOGmin = minimum(minimum.(NANOG))
GATA6max = maximum(maximum.(GATA6))
GATA6min = minimum(minimum.(GATA6))
FGFmax = maximum(maximum.(FGF))
FGFmin = minimum(minimum.(FGF))

nanogmax = alpha/(1+1/(2*K)^(2*coefm))
NANOGmax = nanogmax*Kn
yax = collect(range(0,step=0.01, stop=1))
xax1 = ones(length(yax)).*nth.*NANOGmax
xax2 = ones(length(yax)).*gth.*NANOGmax
xax3 = ones(length(yax)).*nth2.*NANOGmax
xax4 = ones(length(yax)).*gth2.*NANOGmax
histbins = 50

PyPlot.clf()

fig = figure("PDFS",figsize=(20,15)) # Create a new blank figure
subplot(3,4,1)  # used to be normed=True in older versions
x, bins, p=hist(NANOG[1], 3, density=true, color="grey")
y = 0
for item in p
    item.set_height(item.get_height()/sum(x))
end
plot(xax1, yax, linestyle="--", label="EPI th", color="red")
plot(xax2, yax, linestyle="--", label="PrE th", color="blue")
plot(xax3, yax, linestyle="--", label="Upper turn-off", color="black")
plot(xax4, yax, linestyle="--", label="Lower turn-off", color="black")
N = tots[1]
title("N = $N")
ylabel("PDF")
xlabel("NANOG")
ylim(0,1)
xlim(0, NANOGmax)
PyPlot.savefig("PDFS.png")

subplot(3,4,2)
x, bins, p=hist(NANOG[2], histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
plot(xax1, yax, linestyle="--", label="EPI th", color="red")
plot(xax2, yax, linestyle="--", label="PrE th", color="blue")
plot(xax3, yax, linestyle="-.", label="Upper turn-off", color="black")
plot(xax4, yax, linestyle="-.", label="Lower turn-off", color="black")
N = tots[2]
title("N = $N")
xlabel("NANOG")
ylim(0,1)
xlim(0, NANOGmax)

subplot(3,4,3)
x, bins, p=hist(NANOG[3],histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
plot(xax1, yax, linestyle="--", label="EPI th", color="red")
plot(xax2, yax, linestyle="--", label="PrE th", color="blue")
plot(xax3, yax, linestyle="-.", label="Upper turn-off", color="black")
plot(xax4, yax, linestyle="-.", label="Lower turn-off", color="black")
N = tots[3]
title("N = $N")
xlabel("NANOG")
ylim(0,0.5)
xlim(0, NANOGmax)

subplot(3,4,4)
x, bins, p=hist(NANOG[4],histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
plot(xax1, yax, linestyle="--", label="EPI th", color="red")
plot(xax2, yax, linestyle="--", label="PrE th", color="blue")
plot(xax3, yax, linestyle="-.", label="Upper turn-off", color="black")
plot(xax4, yax, linestyle="-.", label="Lower turn-off", color="black")
N = tots[4]
title("N = $N")
xlabel("NANOG")
ylim(0,0.5)
xlim(0, NANOGmax)
PyPlot.legend(bbox_to_anchor=(1.05, 1))

subplot(3,4,5)
x, bins, p=hist(GATA6[1],3, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
ylabel("PDF")
xlabel("GATA6")
ylim(0,1)
xlim(0, GATA6max)

subplot(3,4,6)
x, bins, p=hist(GATA6[2],histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
xlabel("GATA6")
ylim(0,1)
xlim(0, GATA6max)

subplot(3,4,7)
x, bins, p=hist(GATA6[3],histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
xlabel("GATA6")
ylim(0,0.5)
xlim(0, GATA6max)

subplot(3,4,8)
x, bins, p=hist(GATA6[4],histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
xlabel("GATA6")
ylim(0,0.5)
xlim(0, GATA6max)

subplot(3,4,9)
x, bins, p=hist(FGF[1],3, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
ylabel("PDF")
xlabel("FGF")
ylim(0,1)
xlim(0, FGFmax)

subplot(3,4,10)
x, bins, p=hist(FGF[2],histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
xlabel("FGF")
ylim(0,1)
xlim(0, FGFmax)

subplot(3,4,11)
x, bins, p=hist(FGF[3],histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
xlabel("FGF")
ylim(0,0.5)
xlim(0, FGFmax)

subplot(3,4,12)
x, bins, p=hist(FGF[4],histbins, density=true, color="grey")
for item in p
    item.set_height(item.get_height()/sum(x))
end
xlabel("FGF")
ylim(0,0.5)
xlim(0, FGFmax)

plt.tight_layout()

#PyPlot.savefig("PDFS.png")
