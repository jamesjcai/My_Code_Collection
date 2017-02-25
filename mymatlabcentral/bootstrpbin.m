function [m]=bootstrpbin(y)
m = mean(bootstrp(1000,@median,y));