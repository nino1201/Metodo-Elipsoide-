#Numero de variables
#Matriz y vector de condiciones

a=[1 2 3;4 5 6]
f=[10,11]
y=[7,8,9]



#mira para un x si Ax>=b, devuelve I=0 si se cumple, si no devuelve la fila en la que no se cumple
function less(Z,x,b)
	 c=Z*x
	 I=0
	 for i=1:size(x)[1] 
	     if c[i]<b[i]
		I=i
		break
		end
	end
	return I
	end			

#Calcular U
function U(A,b)
	 u1=maximum(abs(A))
	 u2=maximum(abs(b))

	 return max(u1,u2)
	 end

#tiempo para determinar si P es vacio, input V y v con k=1/v 
function t(A,b,n)
	 u=U(A,b)
	 V=(2*n)^n*(n*u)^(n^2)
	 k=n^n*(n*u)^(n^2*(n+1))
	 return ceil(Int64,2*n+1*log(V*k))
	 end


#Calcula el vector x^(k+1) usando x^(k), P(k), gn, n
function evolutionx(x,P,g,n)
	 xn=x+(1/n+1)*P*g/((transpose(g)*P*g)[1])^0.5
	 return xn
	 end


#Calcula la matríz P^(k+1) usando P(k), gn, n	 
function evolutionP(P,g,n)
	 k=n^2/(n^2-1)*(P-2./(n+1)*(P*g*transpose(g)*P)/(transpose(g)*P*g)[1])
	 return k
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
	 #Se calcula el primer "elipsoide", que va a ser una esfera centrada en 0 con radio 2^L/n
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
	       	  print("is empty")
	       	  break
	       end
	 return(x)
	 end
	 end

function SuperM(A,b,c)
	 l=size(A)[1]
	 j=size(A)[2]
	 SM=[A zeros(l,l);zeros(j,j) transpose(A);-transpose(c) transpose(b);transpose(c) -transpose(b);zeros(l,j) eye(l)]
	 return SM
	 end

function SuperVector(b,c)
	 Sv=[b;c;-c;0;0;zeros(size(b)[1])]
	 print(size(Sv))
	 end

function resolver(A,b,c)
	 M=SuperM(A,b,c)
	 h=SuperVector(b,c)
	 print(size(M))
	 end 
	
resolver(a,f,y)