include("MclachlanUpdated.jl")
using Printf


Kop = P-> (a*x+c*x^3)*P+(g*derivative(P,x));
Lop=Dx*(a*X+c*X^3)+(g*Dx^2);
Ldag=-(a*X+c*X^3)*Dx+g*Dx^2;
H = Ldag*Lop;

tfinal=1000.0;
skip = 2;
power = 0:1:4;

numij = 10;
paramtable = [[10.0, 2i+1.0,2j+1.0] for i=0:numij, j=0:numij];
function calcu0s(param_vals::Vector{Float64})
    param_subs = Dict(params::Vector{Num}.=>param_vals)
    Hsubs = substitute(H,param_subs)
    println(param_vals)
    result = calcu0(Hsubs)
    return result[result.>0][1]
end
variationald0 = map(calcu0s,paramtable)
function calcQVMEnergy(H, d0)::Float64
    psi = Hermite(1, -d0, 0, 0)
    psinorm = psi/sqrt(GaussianIntegral(psi^2))

    @assert GaussianIntegral(psinorm^2).val ≈ 1.0

    GaussianIntegral(psinorm*(H*psinorm)).val
end
variationalenergies = [calcQVMEnergy(substitute(H, Dict(params.=>paramtable[i,j])), variationald0[i,j]) for i=1:numij+1, j=1:numij+1]

mclachlanResults=[performMclachlan(i, Lop, params, skip) for i=power];

resultfuncs=transform_result_to_function(mclachlanResults);

solutions = solveEquationsParamsTable(tfinal, H, mclachlanResults, resultfuncs, paramtable, variationald0, flowingTimeSolve, skip);
allenergies = [[variationalenergies];[map(x->x.energy,solutions[:,:,i]) for i=1:length(power)]]

psi11=Hermite(1,-variationald0[1,1],0,0)
psi11norm=psi11/(GaussianIntegral(psi11))
GaussianIntegral(X^4*psi11norm)


stats1=[round.([solutions[1,1,i+1].energy, solutions[1,1,i+1].variance, solutions[1,1,i+1].kurtosis],sigdigits=3) for i=power]
stats1string = ""
for i=1:length(stats1)
    stats1string = stats1string * (i==1 ? "Gaussian" : "$i-Hermite") * " & " * join(stats1[i]," & ") * "\\\\ \\hline \n"
end   
println(stats1string)

#This checks whether the energy of the solutions for the i-Hermite ansatz is lower than the energy of the solutions for the j-Hermite ansatz, or if the i-Hermite ansatz is not NaN and the j-Hermite ansatz is NaN
function isbetterfit(sol1::Float64, sol2::Float64)
    if isnan(sol2)
        return !isnan(sol1)
    elseif isnan(sol1)
        return false
    else
        return sol1 < sol2
    end
end

isbetter = [[[sum(isbetterfit.(allenergies[i],allenergies[j])) for j=1:length(allenergies)];[sum(isnan.(allenergies[i]))]] for i=1:length(allenergies)]
#isgreater=[sum((map(x->x.energy,solutions[:,:,i]).<map(x->x.energy,solutions[:,:,j])))  for i=1:length(power), j=1:length(power)]
isbetterstring=""
for i=1:length(isbetter)
    isbetterstring = isbetterstring * (i==1 ? "QVM Gaussian" : i==2 ? "Gaussian" : "$(i-1)-Hermite") * " & " * join(isbetter[i]," & ") * " \\\\ \\hline \n"
end
println(isbetterstring)

function mean(vec::Vector{T}) where T
    return sum(vec)/length(vec)
end
function variance(vec::Vector{T}) where T
    return sum((vec.-mean(vec)).^2)/(length(vec)-1)
end
function calcAverageDifference(sol1::Matrix{Float64},sol2::Matrix{Float64})
    difference = sol1.-sol2
    nan_indices = Base.:!.(isnan.(difference))
    difference_not_nan = difference[nan_indices]
    weighted_difference = difference_not_nan./sol2[nan_indices]
    return (100*mean(weighted_difference), variance(100 .* weighted_difference))
    #return (100*sum(difference_not_nan./sol2[nan_indices])/length(difference_not_nan), sum(difference_not_nan.<0), length(difference_not_nan))
end
averagedifferenceFromQVM = [calcAverageDifference(allenergies[i+1],variationalenergies) for i=1:length(power)]
averagedifferenceFromGaussian = [calcAverageDifference(allenergies[i+1],allenergies[2]) for i=2:length(power)]
function exportAverageDifference(averagedifference::Vector{Tuple{Float64,Float64}}, names::Vector{String}=String["Gaussian"])
    ret=""
    for i=1:length(averagedifference)
        ret = ret * @sprintf("%.3f",averagedifference[i][1]) * " & "
    end
    ret = ret * " \\\\ \\hline \n"
    for i=1:length(averagedifference)-1
        ret = ret * @sprintf("%.3f",averagedifference[i][2]) * " & "
    end
    ret * @sprintf("%.3f",averagedifference[end][2])
end
println(exportAverageDifference(averagedifferenceFromQVM))

println([sum(isnan.(map(x->x.energy,solutions[:,:,i]))) for i=1:length(power)])
isnan.(map(x->x.energy,solutions[:,:,:]))


## Everything below this line is to export to send to mathematica for plotting Gaussian results
paramtablecsv = [[4.0, 2i+1.0,2j+1.0] for i=0:3, j=0:3];
variationald0csv = map(calcu0s,paramtablecsv);
solutionscsv = solveEquationsParamsTable(tfinal, H, mclachlanResults[1:1], resultfuncs[1:1], paramtablecsv, variationald0csv, flowingTimeSolve, skip);
toexportcsv = map(x->Symbolics.value.([x.ksol[1],-x.dsol[1],x.energy,x.variance,x.kurtosis]),solutionscsv[:,:,1])

weval(W"Export"("Mclachlan Flowing Time Cubic Stats.mx",expr_to_mathematica(toexportcsv)))

map(x->x.energy,solutionscsv[1,:])