#Numero de variables
#Matriz y vector de condiciones

#import JuMP

#mira para un x si Ax>=b, devuelve I=0 si se cumple, si no devuelve la primera fila en la que no se cumple.
function less(Z,x,b)
	 c=Z*x
	 I=0
	 for i=1:size(b)[1] 
	     if c[i]<b[i]
		I=i
		break
		end
	end
	return I
	end

#Calcular U
function U(A,b)
	 An=abs(A)
	 bn=abs(b)
	 u1=maximum(An)
	 u2=maximum(bn)

	 return max(u1,u2)
	 end

#tiempo para determinar si P es vacio, input A, b, n y calcula k=1/v y V.
function t(A,b,n)
	 u=U(A,b)
	 V=log(((2*n)^n)*(n*u)^(n^2))
	 k=log(n^n*(n*u)^((n^2)*(n+1)))
	 return ceil(Int64, 2.*(n+1)*(V+k))
	 end


#Calcula el vector x^(t+1) usando x^(t), P(t), g, n
function evolutionx(x,P,g,n)
	 x = convert(Array{BigFloat}, x)
	 P = convert(Array{BigFloat}, P)
	 xn=x+(1/n+1)*P*g/((transpose(g)*P*g)[1])^0.5
	 return xn
	 end


#Calcula la matríz P^(k+1) usando P(k), gn, n	 
function evolutionP(P,g,n)
	 P = convert(Array{BigFloat}, P)
	 Pn=n^2/(n^2-1)*(P-2./(n+1)*(P*g*transpose(g)*P)/(transpose(g)*P*g)[1])
	 return Pn
	 end


#calcula el primer elipsoide
function e0(A,b,n)
	 P0=n*(n*U(A,b))^(2*n)*eye(n,n)
	 return P0
	 end


#METODO ELIPSOIDE

function ELIPSOIDE(A,b)
	 #dimensión
	 n=size(A)[2]
	 # se halla el número de iteraciones máximas 
	 N=t(A,b,n)
	 # Se calcula el primer elipsoide, que va a ser una esfera centrada en el origen con radio .....
	 P_0=e0(A,b,n)
	 t0=0
	 x=zeros(n)
	 P=P_0
	 while less(A,x,b)!=0
	       i=less(A,x,b)
	       g=transpose(A[i,:])
	       xn=evolutionx(x,P,g,n)
	       Pn=evolutionP(P,g,n) 
	       x=xn
	       P=Pn
	       t0=t0+1
	       if t0==N
	       	  println("is empty \n","\n")
	       	  break
	       end
	 end
	 return(x)
	 end

function SuperM(A,b,c)
	 m=size(A)[1]
	 n=size(A)[2]
	 SM=[A zeros(m,m);zeros(n,n) transpose(A);zeros(n,n) transpose(A);-transpose(c) transpose(b);transpose(c) -transpose(b);zeros(m,n) eye(m)]
	 return SM
	 end

function SuperVector(b,c)
	 Sv=[b;c;-c;0;0;zeros(size(b)[1])]
	 end

function resolver(A,b,c)
	 n=size(A)[2]
	 M=SuperM(A,b,c)
	 h=SuperVector(b,c)
	 e=1/(2(n+1)*((n+1)*U(A,b))^(-(n+1)))*zeros(size(h)[1])
	 h=h-e
	 ELIPSOIDE(M,h)
	 end



#Ejemplo 1: Consiste en verificar si el politopo definido por las desigualdades(x1+x2>=1, x1+x2<=2, x1>=0, x2>=0) es no-vacio
A1=[1. 1.;-1. -1.;1. 0.;0. 1.]
b1=[1.,-2.,0.,0.]

println("Ejemplo 1: ",  ELIPSOIDE(A1,b1) ,"\n")

#Ejemplo 1: con JuMP
#m1=JuMP.Model(solver =ClpSolver())
#JuMP.@variable(m1,x1>=0)
#JuMP.@variable(m1,x2>=0)
#JuMP.@objective(m1,Min,x1+x2)
#JuMP.@constrains(m1,x1+x2>=1)
#JuMP.@constrains(m1,x1+x2<=2)
#print(m1)
#status=solve(m1)
#println("Objective value of Ejemplo1: ", getobjectivevalue(m1))
#println("x1= ",getvalue(x1))
#println("x2= ",getvalue(x2))

#Ejemplo 2: Consiste en verificar si el politopo definido por las desigualdades(x1+2x2<=1, x1+2*x2>=2,x1>=0, x2>=0) es vacio

A2=[-1 -2;1 1;1 0;0 1]
b2=[-1;2;0;0]
print("Ejemplo 2: ")
ELIPSOIDE(A2,b2)
#Ejemplo 2: usando JuMP

#Ejemplo 2: con JuMP
#m2=JuMP.Model(solver =ClpSolver())
#JuMP.@variable(m2,x1>=0)
#JuMP.@variable(m2,x2>=0)
#JuMP.@objective(m2,Min,x1+x2)
#JuMP.@constrains(m2,x1+x2>=2)
#JuMP.@constrains(m2,x1+x2<=1)
#print(m2)
#status2=solve(m2)
#printl(status)

#Ejemplo 3: Consiste en resolver el problema de optimizacion lineal min x1+x2 sujeto a (x1+2*x2<=1, x1+x2>=2,x1>=0, x2>=0)
A3=[1. 0.;0. -1.;1. 0.;0. 1.]
b3=[0.2,-1.,0.,0.]
c=[1.,0.]
print("Ejemplo 3: ",resolver(A1,b1,c))

#Ejemplo 3 Usando JuMP

#m3=JuMP.Model(solver =ClpSolver())
#JuMP.@variable(m3,x1>=0)
#JuMP.@variable(m3,x2>=0)
#JuMP.@objective(m3,Min,x1)
#JuMP.@constrains(m3,x1+x2>=1)
#JuMP.@constrains(m3,x1+x2<=2)
#print(m3)
#status=solve(m3)
#println("Objective value of Ejemplo1: ", getobjectivevalue(m3))
#println("x1= ",getvalue(x1))
#println("x2= ",getvalue(x2))
