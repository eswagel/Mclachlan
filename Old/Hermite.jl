using SpecialFunctions
using SymPy
import SymPy: subs, diff, Latexify.latexify #So I can overload them
SymPy.diff(num::T,sym) where T<:Number = T(0) #Derivative of a number is zero
SymPy.subs(num::Number,rule)=num #Needed for substitution

abstract type HermiteType end
const NumSym=SymPy.Sym#Union{Number,SymPy.Sym} #Eventually get rid of this it causes inefficient multiple dispatch

#This struct models a Hermite basis state
struct Hermite <: HermiteType
    factor::NumSym #Coefficient in front
    a::NumSym #Coefficient of x^2 in the exponent (should be negative)
    b::NumSym #Coefficient of x in the exponent
    order::UInt8 #Power of x in the coefficient
end

Hermite(a::NumSym,b::NumSym,order::Int)=Hermite(1,a,b,order)

#When substituting into a Hermite, substitute into each part 
SymPy.subs(herm::Hermite, args...)::Hermite=Hermite(SymPy.subs(herm.factor, args...), SymPy.subs(herm.a, args...), SymPy.subs(herm.b, args...), herm.order)
#From two arrays, create a Dict for substitution and substitute
SymPy.subs(herm::Hermite, syms1::Union{Array{Sym, 1},Tuple{Sym}}, syms2::Union{Array{Sym, 1},Tuple{Sym}})::Hermite = subs(herm, Dict(syms1[i]=>syms2[i] for i in 1:length(syms1)))
SymPy.subs(herm::Sym, syms1::Union{Vector{Sym},Tuple{Sym}}, syms2::Union{Array{Sym, 1},Tuple{Sym}})::Sym = subs(herm, Dict(syms1[i]=>syms2[i] for i in 1:length(syms1)))

Base.convert(::Type{Sym},herm::Hermite)=herm.factor*x^herm.order*exp(herm.a*x^2+herm.b*x)

#Take the derivative of a single Hermite with respect to some symbol
SymPy.diff(herm::Hermite, var::Sym)::Terms=begin
    dvar = input-> SymPy.diff(input, var)
    result=zero(Terms) #Eventual output
    facdiff=dvar(herm.factor)
    #This if chain essentially does the product rule for factor, a, and b
    if facdiff!=0
        result+=Hermite(facdiff,herm.a,herm.b,herm.order)
    end
    adiff=dvar(herm.a)
    if adiff!=0
        result+=Hermite(herm.factor*adiff,herm.a,herm.b,herm.order+2)
    end
    bdiff=dvar(herm.b)
    if bdiff!=0
        result+=Hermite(herm.factor*bdiff,herm.a,herm.b,herm.order+1)
    end
    return result
end

#Struct that holds an array of Hermites
struct Terms <: HermiteType
    terms::Array{Hermite,1}
end

Terms(herm::Hermite)=Terms([herm])
Terms(terms::Terms)=terms
Base.zero(Terms)=Terms(Hermite[])
Base.zero(::Terms)=Terms(Hermite[]) #Zero is an empty array
Base.zero(::Hermite)=Hermite(0,0,0,0)

#Not working but would be convenient to be able to broadcast over terms
Base.length(terms::Terms)=length(terms.terms)
Broadcast.broadcast(func, terms::Terms) = Terms(Broadcast.broadcast(func, terms.terms))

#Derivative of terms = sum of derivatives
SymPy.diff(terms::Terms, args...)::Terms=begin
    sum=zero(Terms)
    for term in terms.terms
        sum+=SymPy.diff(term, args...)
    end
    sum
end

SymPy.subs(terms::Terms,args...)::Terms=Terms[SymPy.subs.(terms.terms,Ref(args)...)]

Base.convert(::Type{Sym},terms::Terms)=sum(convert.(Ref(Sym),terms.terms))
SymPy.Latexify.latexify(herm::HermiteType,args...;kwargs...)=SymPy.Latexify.latexify(convert(Sym,herm),args...;kwargs...)

#######Group algebra
#Checks if two Hermite basis states can be added together to be one, returns a Bool
addable(herm1::Hermite,herm2::Hermite)::Bool=herm1.a===herm2.a && herm1.b===herm2.b && herm1.order===herm2.order

#Multiplication
Base.:*(herm1::Hermite,herm2::Hermite)::Hermite=Hermite(herm1.factor*herm2.factor,herm1.a+herm2.a,herm1.b+herm2.b,herm1.order+herm2.order)
Base.:*(factor, herm2::Hermite)::Hermite=Hermite(factor*herm2.factor,herm2.a,herm2.b,herm2.order)
Base.:*(herm1::Hermite,factor)::Hermite=factor*herm1

#Exponentiation
Base.:^(herm::Hermite,p)::Hermite=Hermite(herm.factor^p,herm.a*p,herm.b*p,herm.order^p)

#Addition
Base.:+(herm1::Hermite,herm2::Hermite)::Terms=begin
    if herm1.factor==0 && herm2.factor==0
        return zero(Terms)
    elseif addable(herm1,herm2)
        return Terms(Hermite(herm1.factor+herm2.factor, herm1.a, herm1.b, herm1.order))
    elseif herm1.factor==0 || herm1.factor==Sym(0)
        return Terms(herm2)
    elseif herm2.factor==0 || herm2.factor==Sym(0)
        return Terms(herm1)
    else
        return Terms([herm1,herm2])
    end
end
Base.:+(herm::Hermite,terms::Terms)::Terms=begin
    if herm.factor==0
        return Terms(terms.terms)
    end
    #=for i in 1:length(terms.terms)
        if addable(terms.terms[i],herm)
            new_terms=copy(terms.terms)
            new_terms[i]=(terms.terms[i]+herm).terms[1]
            return Terms(new_terms)
        end
    end=#
    return Terms([terms.terms;herm])
end
Base.:+(terms::Terms,herm::Hermite)=herm+terms

Base.:+(a::Terms,b::Terms)=begin
    return Terms([a.terms;b.terms])
end

Base.:-(a::Hermite,b::Terms)::Terms=begin
    sum=Vector{Hermite}(undef,length(b)+1)
    sum[1]=a
    @inbounds @. sum[2:end]=-1*b.terms
    return sum
end
Base.:-(a::Terms,b::Hermite)::Terms=begin
    sum=Vector{Hermite}(undef,length(a)+1)
    @inbounds sum[1:end-1].=a.terms
    @inbounds sum[end]=-b
    return Terms(sum)
end
Base.:-(a::Terms,b::Terms)::Terms=begin
    sum=Vector{Hermite}(undef,length(a)+length(b))
    @inbounds sum[1:length(a)].=a.terms
    @inbounds sum[length(a)+1:end].=-1 .*b.terms
    return Terms(sum)
end
Base.:-(a::HermiteType)=(-1)*a
#Multiplying with Terms is just distribution
Base.:*(a::Terms,b::Terms)::Terms=begin
    new_terms=Vector{Hermite}(undef,length(a)*length(b))
    L=length(b)
    for i=1:length(a)
        for j=1:L
            @inbounds new_terms[(i-1)*L+j]=a.terms[i]*b.terms[j]
        end
    end
    return Terms(new_terms)
end
Base.:*(terms::Terms,herm::Hermite)::Terms=begin
    new_terms=Vector{Hermite}(undef,length(terms))
    @inbounds for i=1:length(terms)
        new_terms[i]=terms.terms[i]*herm
    end
    Terms(new_terms)
end
Base.:*(herm::Hermite,terms::Terms)=terms*herm
Base.:*(fact,terms::Terms)=Terms(fact .* terms.terms)
Base.:*(terms::Terms,fact)=fact*terms

#Integration of a Hermite basis state with a, b, and order.
GaussianIntegral(a, b, order)::Sym=begin
    if b==0
        #If b is 0, we have an explict formula
        return isodd(order) ? zero(Sym) : (-a)^(-(1/2) - order/2)*SpecialFunctions.gamma((1 + order)/2)
    #Otherwise, I've put in some of the lowest orders and could add more or generate a table
    elseif order==0
        result=1/sqrt(-a)
    elseif order==1
        result=-(b/(2*sqrt(-a)*a))
    elseif order==2
        result=-((1 - b^2/(2a))/(2*sqrt(-a)*a))
    elseif order==3
        result=-((3b*(1 - b^2/(6a)))/(4*(-a)^(3/2)*a))
    elseif order==4
        result=(3*(1 - b^2/a + b^4/(12a^2)))/(4*sqrt(-a)*a^2)
    elseif order==5
        result=-((15b*(1 - b^2/(3a) + b^4/(60a^2)))/(8*(-a)^(5/2)*a))
    elseif order==6
        result=-((15*(1 - (3b^2)/(2a) + b^4/(4a^2) - b^6/(120a^3)))/(8*sqrt(-a)*a^3))
    else
        throw(DomainError(order,"Integral for this order Hermite polynomial has not been calculated."))
    end
    return result*sqrt(PI)*exp(-b^2/(4a))
end

GaussianIntegral(herm::Hermite)::Sym=herm.factor*GaussianIntegral(herm.a,herm.b,herm.order)
GaussianIntegral(terms::Terms)::Sym=begin
    ints=Vector{Sym}(undef,length(terms.terms))
    for i in 1:length(terms.terms)
        @inbounds ints[i]=GaussianIntegral(terms.terms[i])
    end
    return sum(ints)
end
struct Operator
    func::Function
end

#Multiplying an operator by a Hermite applies it
Base.:*(op::Operator, herm::Hermite)::Terms=Terms(op.func(herm))
Base.:*(op::Operator, terms::Terms)::Terms=begin
    sum=zero(Terms)
    for i in 1:length(terms.terms)
        sum+=op.func(terms.terms[i])
    end
    return sum
end
#Or the operator can be called directly
(op::Operator)(herm::HermiteType)::Terms=op*herm

#You can also raise an operator to an integer power
Base.:^(op::Operator,expon::Integer)=begin
    result=1
    for j in 1:expon
        result*=op
    end
    return result
end

##Algebra for operators
#For two operators
Base.:*(op1::Operator, op2::Operator)=Operator((herm::Hermite)->op1*(op2*herm))
Base.:+(op1::Operator, op2::Operator)=Operator((herm::Hermite)->(op1*herm+op2*herm))
Base.:-(op1::Operator, op2::Operator)=op1+(-1)*(op2*herm)
Base.:/(op1::Operator, op2::Operator)=op1*(1/(op2*herm))

#For a non-operator and an operator
Base.:*(factor,op::Operator)=Operator((herm::Hermite)->factor*(op*herm))
Base.:+(factor,op::Operator)=Operator((herm::Hermite)->factor+(op*herm))
Base.:-(factor,op::Operator)=Operator((herm::Hermite)->factor-(op*herm))
Base.:/(factor,op::Operator)=Operator((herm::Hermite)->factor/(op*herm))
Base.:-(op::Operator)=(-1)*op

#X raises the order (like multiplying by x)
const x=symbols("x")
function Xfunc(herm::Hermite)::Terms
    Terms(Hermite(herm.factor, herm.a, herm.b, herm.order+1))
end
const X=Operator(Xfunc)
#Derivative with respet to x
function Dxfunc(herm::Hermite)::Terms
    herm.order==0 ? Terms([Hermite(herm.factor*herm.a*2,herm.a,herm.b,1)]) : Terms([Hermite(herm.factor*herm.order,herm.a,herm.b,herm.order-1),Hermite(herm.factor*herm.a*2,herm.a,herm.b,herm.order+1)])
end
const Dx=Operator(Dxfunc)

const t = symbols("t")
