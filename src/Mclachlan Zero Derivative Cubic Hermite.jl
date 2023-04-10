include("MclachlanUpdated.jl")


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

mclachlanResults=[performMclachlan(i, Lop, params, skip) for i=power];

resultfuncs=transform_result_to_function(mclachlanResults);

solutions = solveEquationsParamsTable(tfinal, H, mclachlanResults, resultfuncs, paramtable, variationald0, zeroDerivativesSolve, skip);

size(solutions)

map(x->x.energy,solutions[2,2,:])
map(x->x.variance,solutions[1,1,:])
map(x->x.energy,solutions[:,:,1])-map(x->x.energy,solutions[:,:,2])
solutions[2,2,1].params

@variables P(x)
exactsol = wcall("DSolve",Kop(P)~0,P,x)[1][1][2]
exactnormed = exactsol/wcall("Integrate",exactsol,[x,W`-Infinity`,"Infinity"])
stats1=[round.([solutions[1,1,i+1].energy, solutions[1,1,i+1].variance, solutions[1,1,i+1].kurtosis],sigdigits=3) for i=power]
stats1string = ""
for i=1:length(stats1)
    stats1string = stats1string * (i==1 ? "Gaussian" : "$i-Hermite") * "&" * join(stats1[i],"&") * "\\\\ \\hline \n"
end   
println(stats1string)

isgreater=[sum(map(x->x.energy,solutions[:,:,i]).<map(x->x.energy,solutions[:,:,j])) for i=1:length(power), j=1:length(power)]
isgreaterstring=""
for i=1:length(power)
        isgreaterstring = isgreaterstring * (i==1 ? "Gaussian" : "$i-Hermite") * " & " * join(isgreater[i,:]," & ") * " \\\\ \\hline \n"
end
println(isgreaterstring)

println([sum(isnan.(map(x->x.energy,solutions[:,:,i]))) for i=1:length(power)])


## Everything below this line is to export to send to mathematica for plotting Gaussian results
paramtablecsv = [[4.0, 2i+1.0,2j+1.0] for i=0:3, j=0:3];
variationald0csv = map(calcu0s,paramtablecsv);
solutionscsv = solveEquationsParamsTable(tfinal, H, mclachlanResults[1:1], resultfuncs[1:1], paramtablecsv, variationald0csv, zeroDerivativesSolve, skip);
toexportcsv = map(x->Symbolics.value.([x.ksol[1],-x.dsol[1],x.energy,x.variance,x.kurtosis]),solutionscsv[:,:,1])

weval(W"Export"("Mclachlan Zero Derivative Cubic.mx",expr_to_mathematica(toexportcsv)))

map(x->x.energy,solutionscsv[:,:])