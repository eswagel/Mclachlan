using Distributed
using SpecialFunctions
using SymPy
import SymPy: subs, diff, Latexify.latexify #So I can overload them
SymPy.diff(num::T,sym) where T<:Number = T(0) #Derivative of a number is zero
SymPy.subs(num::Number,rule)=num #Needed for substitution

abstract type HermiteType end

#This struct models a Hermite basis state
struct Hermite <: HermiteType
    factor::SymPy.Sym #Coefficient in front
    a::SymPy.Sym #Coefficient of x^2 in the exponent (should be negative)
    b::SymPy.Sym #Coefficient of x in the exponent
    order::UInt8 #Power of x in the coefficient
end

Hermite(a::SymPy.Sym,b::SymPy.Sym,order::Int)=Hermite(1,a,b,order)

#When substituting into a Hermite, substitute into each part 
subs(herm::Hermite, args...)::Hermite=Hermite(SymPy.subs(herm.factor, args...), SymPy.subs(herm.a, args...), SymPy.subs(herm.b, args...), herm.order)
#From two arrays, create a Dict for substitution and substitute
subs(herm::Hermite, syms1::Union{Array{SymPy.Sym, 1},Tuple{SymPy.Sym}}, syms2::Union{Array{SymPy.Sym, 1},Tuple{SymPy.Sym}})::Hermite = subs(herm, Dict(syms1[i]=>syms2[i] for i in 1:length(syms1)))
subs(herm::SymPy.Sym, syms1::Union{Vector{SymPy.Sym},Tuple{SymPy.Sym}}, syms2::Union{Array{SymPy.Sym, 1},Tuple{SymPy.Sym}})::SymPy.Sym = subs(herm, Dict(syms1[i]=>syms2[i] for i in 1:length(syms1)))

convert(::Type{SymPy.Sym},herm::Hermite)=herm.factor*x^herm.order*exp(herm.a*x^2+herm.b*x)

#Take the derivative of a single Hermite with respect to some symbol
SymPy.diff(herm::Hermite, var::SymPy.Sym)::Terms=begin
    dvar = input-> SymPy.diff(input, var)
    result=zero(Terms) #Eventual output
    facdiff=dvar(herm.factor)
    #This if chain essentially does the product rule for factor, a, and b
    if facdiff!=0
        @inbounds result+=Hermite(facdiff,herm.a,herm.b,herm.order)
    end
    adiff=dvar(herm.a)
    if adiff!=0
        @inbounds result+=Hermite(herm.factor*adiff,herm.a,herm.b,herm.order+2)
    end
    bdiff=dvar(herm.b)
    if bdiff!=0
        @inbounds result+=Hermite(herm.factor*bdiff,herm.a,herm.b,herm.order+1)
    end
    return result
end

#Struct that holds an array of Hermites
struct Terms <: HermiteType
    terms::Array{Hermite,1}
end

Terms(herm::Hermite)=Terms([herm])
Terms(terms::Terms)=terms
Base.zero(Terms)=Terms([]) #Zero is an empty array

#Not working but would be convenient to be able to broadcast over terms
Base.length(terms::Terms)=length(terms.terms)
Broadcast.broadcast(func, terms::Terms) = Terms(Broadcast.broadcast(func, terms.terms))

#Derivative of terms = sum of derivatives
SymPy.diff(terms::Terms, args...)::Terms=begin
    sum=zero(Terms)
    @sync @distributed for term in terms.terms
        sum+=SymPy.diff(term, args...)
    end
    sum
end

subs(terms::Terms,args...)::Terms=begin
    for i=1:length(terms.terms)
        terms.terms[i]=subs(terms.terms[i],args...)
    end
    return terms
end

convert(::Type{SymPy.Sym},terms::Terms)=sum(convert.(Ref(SymPy.Sym),terms.terms))
SymPy.Latexify.latexify(herm::HermiteType,args...;kwargs...)=SymPy.Latexify.latexify(convert(SymPy.Sym,herm),args...;kwargs...)
Sym(herm::HermiteType)=convert(SymPy.Sym,herm)

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
    elseif herm1.factor==0 || herm1.factor==SymPy.Sym(0)
        return Terms(herm2)
    elseif herm2.factor==0 || herm2.factor==SymPy.Sym(0)
        return Terms(herm1)
    else
        return Terms([herm1,herm2])
    end
end
Base.:+(herm::Hermite,terms::Terms)::Terms=begin
    if herm.factor==0
        return terms
    end
    #=for i in 1:length(terms.terms)
        if addable(terms.terms[i],herm)
            new_terms=copy(terms.terms)
            @inbounds new_terms[i]=(terms.terms[i]+herm).terms[1]
            return Terms(new_terms)
        end
    end=#
    push!(terms.terms,herm);
    terms
end
Base.:+(terms::Terms,herm::Hermite)=herm+terms

Base.:+(a::Terms,b::Terms)=begin
    sum=a
    @sync @distributed for i in 1:length(b.terms)
        @inbounds sum+=b.terms[i]
    end
    sum
end

#Subtraction means just multiplying by -1
Base.:-(a::HermiteType,b::HermiteType)=a+(-1)*b
Base.:-(a::HermiteType)=(-1)*a

#Multiplying with Terms is just distribution
Base.:*(a::Terms,b::Terms)::Terms=begin
    new_terms=zero(Terms)
    @sync @distributed for term in a.terms
        @inbounds new_terms+=term*b
    end
    new_terms
end
Base.:*(terms::Terms,herm::Hermite)::Terms=begin
    @sync @distributed for i=1:length(terms.terms)
        @inbounds terms.terms[i]=terms.terms[i]*herm
    end
    terms
end
Base.:*(herm::Hermite,terms::Terms)=terms*herm
Base.:*(fact,terms::Terms)=Terms(fact .* terms.terms)
Base.:*(terms::Terms,fact)=fact*terms

Base.:/(herm::HermiteType,fact)=herm*(1/fact)

#Not used but can be useful to define
Base.sum(terms::Terms)=terms
Base.sum(terms::Array{HermiteType,1})=Terms(terms)

#Integration of a Hermite basis state with a, b, and order.
GaussianIntegral(a, b, order)::SymPy.Sym=begin
    if b==0
        #If b is 0, we have an explict formula
        return isodd(order) ? 0 : (-a)^(-(0.5) - order/2)*SpecialFunctions.gamma((order+1)/2)
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

GaussianIntegral(herm::Hermite)::SymPy.Sym=herm.factor*GaussianIntegral(herm.a,herm.b,herm.order)
GaussianIntegral(terms::Terms)::SymPy.Sym=begin
    ints=SymPy.Sym(0)
    @sync @distributed for i in 1:length(terms.terms)
        ints+=GaussianIntegral(terms.terms[i])
    end
    return ints
end

struct Operator
    func::Function
end

#Multiplying an operator by a Hermite applies it
Base.:*(op::Operator, herm::Hermite)::Terms=op.func(herm)
Base.:*(op::Operator, terms::Terms)::Terms=begin
    sum=zero(Terms)
    @sync @distributed for i in 1:length(terms.terms)
        @inbounds sum+=op.func(terms.terms[i])
    end
    return sum
end
#Or the operator can be called directly
(op::Operator)(herm::HermiteType)=op*herm

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
const X=Operator((herm::Hermite)->Terms(Hermite(herm.factor, herm.a, herm.b, herm.order+1)))
#Derivative with respet to x
const Dx=Operator((herm::Hermite)->herm.order==0 ? Terms([Hermite(herm.factor*herm.a*2,herm.a,herm.b,1)]) : Terms([Hermite(herm.factor*herm.order,herm.a,herm.b,herm.order-1),Hermite(herm.factor*herm.a*2,herm.a,herm.b,herm.order+1)]))
const Dt=Operator((herm::Hermite)->Terms(SymPy.diff(herm,t)))

const t = symbols("t")
