function value = roundn(x, n)

 % Round to nearest with 10^n precision
 value = round(x/(10^n))*(10^n);
 
end
