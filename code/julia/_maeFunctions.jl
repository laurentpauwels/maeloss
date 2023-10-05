"""
# _maeFunctions.jl 

This jl file contains all the functions required to reproduce all the results in the Jupyter Notebooks. 
- This file is NOT necessary to run the Jupyter Notebook. 
- This is a collection of the main functions for convenience.

## Structure 
There are three main sections in this file covering the construction of Figure 1(a) and 1(b), Simulation results and 
the empirical resutls. Packages required for each section will be loaded at the beginning of each section. 
For the overall requirements of all the packages, see requirement.txt. 


"""

######################################################################### Section 1. Figure 1(a) and 1(b) #################################################################################################################################

using QuadGK, Distributions

##############################################################################################################################################################################################################################################
"""
## Usage
skewNormalMSE(a::Float64, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Float64

## Description:

Calcualte the expected value of `z = \\alpha y_1 + (1-\\alpha) y_2` when `y_i` follows an independent skewed normal distribution. Note that in this function, the restriction that the sum of the two weights equals to 1 is imposed. 

The skewed normal distribution is presumed to the following stochastic representationn 

```{math}
y_i = \\xi_i + \\lambda_i \\tau_i + u_i
```

where `\\tau_i` and `u_i` follow a half normal and a normal respectively. 

## Input:
    a: Float 64. α weight.
    lambda: Array{Float64, 1}(2,1). Skewed parameters for the two independent skewed normal random variables. 
    sigma: Array{Float64, 1}(2,1). The variance paramters in `u_i`. 

## Output:
    Float64. The expected value of `z^2`. 

"""
function skewNormalMSE(a::Float64, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Float64
    delta = a*lambda[1] + (1-a)*lambda[2]
    tsigma = a^2*sigma[1] + (1-a)^2*sigma[2]
    return delta^2*(1-2/pi) + tsigma
end

##############################################################################################################################################################################################################################################

"""
## Usage 
    SkewNormalExpectedAbsolute(a::Float64, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Float64

## Description

Calculate the expected value of `|z|` where `z = \\alpha y_1 + (1-\\alpha) y_2` with `y_i` following an independent skewed normal distribution. Note that in this function, the restriction that the sum of the two weights equals to 1 is imposed.  
The skewed normal distribution is presumed to the following stochastic representationn 

```{math}
y_i = \\xi_i + \\lambda_i \\tau_i + u_i
```

where `\\tau_i` and `u_i` follow a half normal and a normal respectively. 

## Input:
    a: Float 64. α weight.
    lambda: Array{Float64, 1}(2,1). Skewed parameters for the two independent skewed normal random variables. 
    sigma: Array{Float64, 1}(2,1). The variance paramters in `u_i`. 

## Output:
    Float64. The expected value of `|z|`. 

"""
function SkewNormalExpectedAbsolute(a::Float64, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Float64
    delta = a*lambda[1] + (1-a)*lambda[2]
    tsigma = a^2*sigma[1] + (1-a)^2*sigma[2]
    densityform = Distributions.Normal(-delta*(2/pi)^(0.5), tsigma)
    cdfform = Distributions.Normal(0, 1-(delta^2)/tsigma)
    density(x) = abs(x)*Distributions.pdf(densityform, x)*Distributions.cdf(cdfform, (x+delta*(2/pi)^(0.5))*delta/tsigma)
    integrate = quadgk(density, -Inf, Inf)
    return integrate[1]
end


##############################################################################################################################################################################################################################################
#
"""
## Usage 

    skewNormalMAE(a::Array{Float64,1}, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Array{Float64}

## Description
    
    A wrapper function of *SkewNormalExpectedAbsolute* allowing a vector of α as input. *SkewNormalExpectedAbsolute* applies to each element in the α vector with the same ξ, λ and σ inputs.  

## Input:
    a: Array{Float64, 1}. Vector of α weights.
    lambda: Array{Float64, 1}(2,1). Skewed parameters for the two independent skewed normal random variables. 
    sigma: Array{Float64, 1}(2,1). The variance paramters in `u_i`. 

## Output
    An Array{Float64,1} the same size as a containing the expected values of `|z|` for each value in a. 
"""
function skewNormalMAE(a::Array{Float64,1}, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Array{Float64}
    temp(s) = SkewNormalExpectedAbsolute(s, lambda, sigma)
    return temp.(a)
end

##############################################################################################################################################################################################################################################

"""
## Usage
skewNormalMSE(a::Array{Float64,1}, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Float64

## Description:

Calcualte the expected value of `z = \\alpha y_1 + (1-\\alpha) y_2` when `y_i` follows an independent skewed normal distribution. Note that in this dispatch, there is no restriction on the sum of the weight. The first argument requires a `2 \\times 1` array containing the values of `\\alpha_1` and `\\alpha_2`.  

The skewed normal distribution is presumed to the following stochastic representationn 

```{math}
y_i = \\xi_i + \\lambda_i \\tau_i + u_i
```

where `\\tau_i` and `u_i` follow a half normal and a normal respectively. 

## Input:
    a: Array{Float64}(2,1). `\\alpha_1` and `\\alpha_2` weights.
    lambda: Array{Float64, 1}(2,1). Skewed parameters for the two independent skewed normal random variables. 
    sigma: Array{Float64, 1}(2,1). The variance paramters in `u_i`. 

## Output:
    Float64. The expected value of `z^2`. 

"""
function skewNormalMSE(a::Array{Float64,1}, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Float64
    delta = transpose(a)*lambda
    tsigma = transpose(a.^2)*sigma
    return delta^2*(1-2/pi) + tsigma
end


##############################################################################################################################################################################################################################################
#
"""
## Usage 
    SkewNormalExpectedAbsolute(a::Array{Float64,1}, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Float64

## Description

Calculate the expected value of `|z|` where `z = \\alpha y_1 + (1-\\alpha) y_2` with `y_i` following an independent skewed normal distribution. Note that in this dispatch, there is no restriction on the sum of the two weights and the first argument requires a `2 \\times 1` array containing the values of both `\\alpha_1` and `\\alpha_2`. 

The skewed normal distribution is presumed to the following stochastic representationn 

```{math}
y_i = \\xi_i + \\lambda_i \\tau_i + u_i
```

where `\\tau_i` and `u_i` follow a half normal and a normal respectively. 

## Input:
    a: Array{Float64,1}. `\\alpha_1` and `\\alpha_2` weights.
    lambda: Array{Float64, 1}(2,1). Skewed parameters for the two independent skewed normal random variables. 
    sigma: Array{Float64, 1}(2,1). The variance paramters in `u_i`. 

## Output:
    Float64. The expected value of `|z|`. 

"""
function SkewNormalExpectedAbsolute(a::Array{Float64,1}, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Float64
    delta = transpose(a)*lambda
    tsigma = transpose(a.^2)*sigma
    densityform = Distributions.Normal(-delta*(2/pi)^(0.5), tsigma)
    cdfform = Distributions.Normal(0, 1-(delta^2)/tsigma)
    density(x) = abs(x)*Distributions.pdf(densityform, x)*Distributions.cdf(cdfform, (x+delta*(2/pi)^(0.5))*delta/tsigma)
    integrate = quadgk(density, -Inf, Inf)
    return integrate[1]
end

##############################################################################################################################################################################################################################################

"""
## Usage 

    skewNormalMAE(a::Array{Float64,2}, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Array{Float64}

## Description
    
    A wrapper function of *SkewNormalExpectedAbsolute* allowing a vector of α as input without restriction on the sum of the weights. *SkewNormalExpectedAbsolute* applies to each element in the α vector with the same ξ, λ and σ inputs.  

## Input:
    a: Array{Float64, 2}(N,2). Matrix of α weights where .
    lambda: Array{Float64, 1}(2,1). Skewed parameters for the two independent skewed normal random variables. 
    sigma: Array{Float64, 1}(2,1). The variance paramters in `u_i`. 

## Output
    An Array{Float64,1} the same size as a containing the expected values of `|z|` for each value in a. 
"""
function skewNormalMAE(a::Array{Float64,2}, lambda::Array{Float64,1}, sigma::Array{Float64,1})::Array{Float64}
    temp(s) = SkewNormalExpectedAbsolute(s, lambda, sigma)
    return [temp(convert(Array{Float64,1},sa)) for sa in eachrow(a)]
end

##############################################################################################################################################################################################################################################


########################################################################### Section 2: Simluations ###########################################################################################################################################

"""
## Usage
    MSEOptimal(Omega::AbstractArray{Float64, 2}; inverse=false)::Array{Float64,1} 

## Description
    Calculate the MSE optimal weight given the variance-covariance matrix of the forecast errors, `\\Omega`. If inverse is true then it will calculate the inverse of `\\Omega` directly using the the LinearAlgebra package. Otherwise, the inverse will be obtained as the solution to the simultaneous equation `\\mathbf{A} \\mathbf{x} = \\mathbf{I}` using the operator "/".

## Inputs
    Omega: AbstractArray{Float64,2}. Variance-covariance matrix of the forecast errors. 
    inverse: Boolean. If true, it uses the inverse function from LinearAlgebra to compute the inverse of Omega. Otherwise, it uses the solution to a simultaneous equation (see description above). 

## Output
    a: Array{Float64,1). The vector of optimal MSE weight. 

"""
function MSEOptimal(Omega::AbstractArray{Float64, 2}; inverse=false)::Array{Float64,1} 
    p = size(Omega,1)
    i = ones(p)
    if inverse == false
    	a = Omega\ones(p,1)/(ones(p,1)'/Omega*ones(p,1));
    else
	invOmega = inv(Omega)
	H = transpose(i)*invOmega*i
	a = invOmega*i/H
    end
    return reshape(a,p) 
end

"""
"""
function MAEOptimal(Z::AbstractArray{Float64,2})::Array{Float64,1}
    T,p = size(Z)
    Ω = transpose(Z)*Z/T
    a0 = MSEOptimal(Ω)
    mae_optim = optimize(a->maeloss(a,Ω), a0, BFGS())
    a = minimizer(mae_optim)
    return reshape(a,p)
end


"""
## Usage 
    maeloss(a_init::Array{Float64,1},nu_it::Array{Float64,2})::Float64

## Description

    Calculate MAE Lost given the weight vector and the forecast errors. The weights are assumed to be summed to 1. 

## Inputs
    
    a_init: Array{Float64, 1}. `(p-1)\\times 1` vector containing the weights. The weight vector is assumed to be summed to 1 so the `p^{th}` element will be calculated as 1 substract the sum of a_init. 

    nu_it:: Array{Float64, 2}. `T\\times p` matrix containing `T` forecast errors from `p` different models. 

## Output

    mae: Float64. The MAE from the `p` models using the weight vector. 

"""
function maeloss(a_init::Array{Float64,1},nu_it::Array{Float64,2})::Float64
    (T,p) = size(nu_it)
    a = zeros(p,1)
    a[1:p-1,1] = a_init
    a[p,1] = 1 - sum(a_init)
    nu_c = nu_it*a
    mae = sum(abs.(nu_c))/T
    return mae
end


"""
"""
function SimulateSkewedNormal(ξ::AbstractArray{Float64,1}, Lambda::AbstractArray{Float64}, Σ::Abstract{Float64,2}, period::Int)::Array{Float64,2}
    p = size(Σ,1)
    tau = abs.(rand(Normal(0,1), (p,period)))
    U = rand(MvNormal(Σ), period)
    if size(Lambda,2) == 1
	Λ = Diagonal(Lambda)
    else
	Λ = Lambda
    y = kronecker(ones(1,period), ξ) + Λ*tau + U 
    return y
end

"""
"""
function simulations(atheory,Sigma,Lambda,R,N)
    p = size(Sigma,1)
    n = length(N)
    asn_MSE = Array{Float64,3}(undef, (p,R,n)) 
    asn_MAE = Array{Float64,3}(undef, (p,R,n))  
    at3_MSE = Array{Float64,3}(undef, (p,R,n)) 
    at3_MAE = Array{Float64,3}(undef, (p,R,n)) 
    for (i, Ni) in enumerate(N) 
        for r = 1:R
        # Skew Random Forecast Errors
            Xi = -sqrt(2/pi)*Lambda*ones(p,1)
	    Z = SimulateSlewedNormal(reshape(Xi,p), Lambda Sigma, N) 
        
        #MSE weights
            Omegahat_SN = (Z'*Z)/Ni
            asn_MSE[:,r,i] = Omegahat_SN\ones(p,1)/(ones(p,1)'/Omegahat_SN*ones(p,1))
            #display(asn_MSE
        #MAE weights
            asn_init = Omegahat_SN\ones(p,1)/(ones(p,1)'/Omegahat_SN*ones(p,1))
            sln_sn = optimize(asn_mae->maeloss(asn_mae, Z), asn_init[1:p-1,1],BFGS()) 
            asn_mae_raw = Optim.minimizer(sln_sn)
            asn_MAE[:,r,i] = [asn_mae_raw;1-sum(asn_mae_raw)]
        
        ## t_3 forecast errors
            Omega = Sigma + (1-2/pi)*(Lambda*Lambda')
            E = rand(MvTDist(df,zeros(p),Omega),20)'
        
        #MSE weights
            Omegahat_t3 = (E'*E)/Ni;
            at3_MSE[:,r,i] = Omegahat_t3\ones(p,1)/(ones(p,1)'/Omegahat_t3*ones(p,1))
        
        #MAE weights
            at3_init = Omegahat_t3\ones(p,1)/(ones(p,1)'/Omegahat_t3*ones(p,1))
            sln_t3 = optimize(at3_mae->maeloss(at3_mae, Z), at3_init[1:p-1,1],BFGS()) 
            at3_mae_raw = Optim.minimizer(sln_t3)
            at3_MAE[:,r,i] = [at3_mae_raw;1-sum(at3_mae_raw)];
            #display(at3_MAE)
        end
    end
    display("Simulations are done")
    return asn_MSE, asn_MAE, at3_MSE, at3_MAE
end


##################################################################### Section 3: Empirical Illustration ##########################################################################################################################################


"""
## Usage 
    expandingwindow(y::Array{Float64,1},Z::Array{Float64,2}; portion::Float64=0.5)::Tuple{AbstractArray{Float64,2}, AbstractArray{Float64,2}}

## Descriptions

    Function to calculate the optimal weight based on recursive windows. 

## Inputs

    y: Array{Float64,1}. A `T\\times 1` vector containing the actual values 
    Z: Array{Float64,2}. `T\\times p` matrix containing the forecast values
    
"""
function expandingwindow(y::AbstractArray{Float64},Z::AbstractArray{Float64}; portion::Float64=0.5)::Tuple{AbstractArray{Float64}, AbstractArray{Float64}}
    (T,p) = size(Z)
    T2 = Int.(collect(range(floor(T*portion), stop=T, step=1)))
    a_MSE = zeros(length(T2),p)
    a_MAE = zeros(length(T2),p)
    allerrors = kron(ones(1,p), y) - Z
    for i = 1:length(T2)
            t = T2[i]
	    errors = allerrors[1:t,:]
            #MSE weights 
            Omegahat = cov(errors);
	    invOmegahat = inv(Omegahat)
	    #a_MSE[i,:] = reshape(invOmegahat*ones(p,1)*inv(ones(1,p)*invOmegahat*ones(p,1)), (1,p)) #using inv produces slightly different results but still qualitative consistent. 
	    a_MSE[i,:] = reshape(Omegahat\ones(p,1)/(ones(p,1)'/Omegahat*ones(p,1)), (1,p))
	    #MSE weights
	    a_init = reshape(a_MSE[i,:], p)
            sln = optimize(a_mae->maeloss(a_mae, errors), a_init[1:p-1],BFGS()) 
            a_mae_raw = Optim.minimizer(sln)
            a_MAE[i,:] = [a_mae_raw;1-sum(a_mae_raw)]   
    end
    return a_MSE, a_MAE
end

##############################################################################################################################################################################################################################################

"""

## Usage

    getHighMomentStats(y::AbstractArray{Float64},Z::AbstractArray{Float64})::Tuple{AbstractArray{Float64,1}, AbstractArray{Float64,1}, AbstractArray{Float64,1}}

## Descriptions
    
    Compute skewness, kurtosis and JB test statistics for a set of forecast errors. Both skewness and kurtosis are defined as standardised central third and fourth moments. 

## Inputs
   
   y::AbstractArray{Float64}. A `T\\times 1` vector containing the actual values of the target variable. 
   Z::AbstractArray{Float64}. A `T\\times p` matrix containing the forecasts of `y` from `p` different models. 

## Output
    A tuple of the form
	Tuple{AbstractArray{Float64,1}, AbstractArray{Float64,1}, AbstractArray{Float64,1}}
    which elements are the Kurtosis, Skewness and the p-value of the JB test statistics, respectively.  

"""
function getHighMomentStats(y::AbstractArray{Float64},Z::AbstractArray{Float64})::Tuple{AbstractArray{Float64,1}, AbstractArray{Float64,1}, AbstractArray{Float64,1}}
    (T,p) = size(Z)
    ferrors = kron(ones(1,p),y) - Z
    errors_kurt = Array{Float64,1}(undef, p) 
    errors_skew = Array{Float64,1}(undef, p) 
    jbpval = Array{Float64,1}(undef, p) 
    for i = 1:p
        errors_kurt[i] = 3+kurtosis(ferrors[:,i])
        errors_skew[i] = skewness(ferrors[:,i])
        jbpval[i] = pvalue(JarqueBeraTest(ferrors[:,i]))
    end
    return errors_kurt, errors_skew, jbpval
end

##############################################################################################################################################################################################################################################

"""
## Usage 
    getExpandWindowsForecastCriteria(y::AbstractArray{Float64},Z::AbstractArray{Float64},a_MSE::AbstractArray{Float64,2},a_MAE::AbstractArray{Float64,2}; portion::Float64=0.5)::Tuple{Array{Float64,1}, Array{Float64,1}}

## Description
    Calculate the Mean Squared Errors and Mean Absolute Errors from a set of weights given by the function *expandingwindow*. The weights obtained from 1 to `t` observations will be used to combined the `(t+1)^th` forecasts to predict the `t+1` observation. 

## Inputs
    y::AbstractArray{Float54}. A `T\\times 1` vector contianing the actual values of the target variable. 
    Z::AbstractArray{Float64}. A `T\\times p` matrix containing `T` forecasts from `p` models. 
    a_MSE::AbstractArray{Float64, 2}. A `n \\times p` matrix containing `n` sets of MSE optimal weights for the `p` models where `n` equals portion (see below). This is most likely to be the output from *expandwindow* (see expandwindow). 
    a_MAE::AbstractArray{Float64,2}. A `n \\times p` matrix containing `n` sets of MAE optimal weights for the `p` models where `n` equals portion (see below). This is most likely to be the output from *expandwindow* (see expandwindow).  
    portion::Float64. The portion of sample to obtain the first set of weight. 

"""
function getExpandWindowsForecastCriteria(y::AbstractArray{Float64},Z::AbstractArray{Float64},a_MSE::AbstractArray{Float64,2},a_MAE::AbstractArray{Float64,2}; portion::Float64=0.5)::Tuple{Array{Float64}, Array{Float64}}
    (T,p) = size(Z)
    T2 = Int(floor(T*portion))
    errors_wmse = y[T2+1:end,:] - sum(a_MSE[1:end-1,:].*Z[T2+1:end,:], dims=2)
    errors_wmae = y[T2+1:end,:] - sum(a_MAE[1:end-1,:].*Z[T2+1:end,:], dims=2)
    Npred = size(errors_wmse,1);
    MSFE = [sum(errors_wmse.^2,dims=1)/Npred;sum(errors_wmae.^2,dims=1)/Npred]
    MAFE = [sum(abs.(errors_wmse),dims=1)/Npred;sum(abs.(errors_wmae),dims=1)/Npred]
    return MSFE, MAFE
end
