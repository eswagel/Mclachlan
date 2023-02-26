const sympy_parse_mathematica = SymPy.PyCall.pyimport("sympy.parsing.mathematica")
using MathLink
import MathLink: W2Mstr

W2Mstr(x::MathLink.WReal)=begin
    return SubString(string(x),findfirst('\`',string(x)))
end
function W2Mstr_PLUS(x::Union{Array,Tuple})
    #println("W2Mstr_PLUS:",x)
    Str="("
    for j in 1:length(x)
        if j>1
            Str*=" + "
        end
        Str*=W2Mstr(x[j])
    end
    Str*=")"
end
    
function W2Mstr_TIMES(x::Union{Array,Tuple})
    #println("W2Mstr_TIMES:",x)
    Str="("
    for j in 1:length(x)
        if j>1
            Str*="*"
        end
        Str*=W2Mstr(x[j])
    end
    Str*=")"
end
    
    
function W2Mstr_POWER(x::Union{Array,Tuple})
    #println("W2Mstr_POWER:",x)
    if length(x) != 2
        error("Power takes two arguments")
    end
    Str="("*W2Mstr(x[1])*"^"*W2Mstr(x[2])*")"
end


    
function W2Mstr_COMPLEX(x::Union{Tuple,Array})
    #println("W2Mstr_COMPLEX:",x)
    if length(x) != 2
        error("Complex takes two arguments")
    end
    if x[1] == 0
        ###Imaginary
        Str="("*W2Mstr(x[2])*"*I)"
    elseif x[2] == 0
        ### Real
        ###Complex
        Str=W2Mstr(x[1])
    else
        ###Complex
        Str="("*W2Mstr(x[1])*"+"*W2Mstr(x[2])*"*I)"
    end
end
function W2Mstr(x::MathLink.WExpr)
    #println("W2Mstr::",x.head.name)
    if x.head.name == "Plus"
        Str = W2Mstr_PLUS(x.args)
    elseif x.head.name == "Times"
        Str = W2Mstr_TIMES(x.args)
    elseif x.head.name == "Power"
        Str = W2Mstr_POWER(x.args)
    elseif x.head.name == "Complex"
        Str = W2Mstr_COMPLEX(x.args)
    elseif x.head.name == "Rational"
        Str = string(x.args[1]/x.args[2])
    else
        Str=x.head.name*"["
        for j in 1:length(x.args)
            if j>1
                Str*=","
            end
            Str*=W2Mstr(x.args[j])
        end
        Str*="]"
    end
    return Str
end

DSolveMathematica(expr::AbstractVector{SymPy.Sym}, vars::AbstractVector{SymPy.Sym}; ics=Dict())::Vector{SymPy.Sym}=DSolveMathematica(((expr)->SymPy.sympy.printing.mathematica.mathematica_code(expr.lhs)*"=="*SymPy.sympy.printing.mathematica.mathematica_code(expr.rhs)).(expr), vars; ics=ics)
DSolveMathematica(expr::SymPy.Sym, vars::AbstractVector{SymPy.Sym}; ics=Dict())::Vector{SymPy.Sym}=DSolveMathematica(SymPy.sympy.printing.mathematica.mathematica_code(expr), vars; ics=ics)
function DSolveMathematica(expr::Vector{String}, vars::AbstractVector{SymPy.Sym}; ics=Dict())::Vector{SymPy.Sym}
    mathVars::String="{" * join(SymPy.sympy.printing.mathematica.mathematica_code.(vars),",") * "}"
    exprWithIcs::String="{"*join(expr,",")
    for (key, value) in ics
        exprWithIcs*="," * SymPy.sympy.printing.mathematica.mathematica_code(key) * "==" * SymPy.sympy.printing.mathematica.mathematica_code(value)
    end
    exprWithIcs*="}"
    mathExpr::MathLink.WExpr=W`N`(W`DSolveValue`(W`ReleaseHold`(W`$exprWithIcs`),W`$mathVars`,W`t`))
    dsolved::MathLink.WExpr=weval(mathExpr)
    len::Integer=weval(W`Length`(dsolved))
    println([W2Mstr(weval(W`Part`(dsolved,i))::MathLink.WExpr) for i=1:len])
    [SymPy.sympify(W2Mstr(weval(W`Part`(dsolved,i))::MathLink.WExpr)::String)::Sym for i=1:len]
end

