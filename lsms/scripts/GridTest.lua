print("Testing the Radial Grid for LSMS_3")
x0=0.1
h=0.1
jws=501
print("Generating a grid with x0 = ",x0,", h=",h," and ",jws," points")
r=RadialGrid.new()
r:generate(x0,h,jws,jws,jws)
for i=0, r:size()-1 do
  print(i," ",r:x(i)," ",r[i])
end