function r = pdf_truncnorm(x,mu,sigma) 

r = normpdf(x,mu,sigma) ./ (1-normcdf(-1,mu,sigma));
