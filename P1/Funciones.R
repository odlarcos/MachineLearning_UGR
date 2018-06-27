# Ejercicio 1

gradientDescent <- function(w,n,maxiters,epsilon,funcion,derivada){
  
  X <- w
  Y <- funcion(X)
  
  error <- (1/N)*( norm((X%*%w)-Y,"F")^2 )
  i=1
  while( (i < maxiters) & (epsilon < error) ){
    w = w - n * derivada(w)
    
    error <- (1/N)*( norm((X%*%w)-Y,"F")^2 )
    i=i+1
  }
  w
}

# Ejercicio 2

E <- function(u,v){
  (((u^3)*exp(v-2))-(4*(v^3)*exp(-u)))^2
}

DerivadaE <- function(w){
  u = w[1]
  v = w[2]
  
  du = 2 * ((3 * u^2 * exp(v - 2) + 4 * (v^3) * exp(-u)) * (((u^3) * exp(v - 2)) - (4 * (v^3) * exp(-u))))
  dv = 2 * (((u^3) * exp(v - 2) - 4 * (3 * v^2) * exp(-u)) * (((u^3) * exp(v - 2)) - (4 * (v^3) * exp(-u))))
  
  c(du,dv)
}

gradientDescent <- function(w,n,maxiters,funcion,derivada){
  
  control <- 1000000
  
  i=1
  fin = FALSE
  while(i < maxiters & !fin){
    w = w - n * derivada(w)
    
    if(i == control){
      sprintf("Iteracion: %i \n", i)
      control <- (control+control)
    }
    
    if(funcion(w[1],w[2]) < 10^(-14)){
      result = c(i, w)
      fin = TRUE
    }
    i=i+1
  }
  names(result) <- c("Iteración","Coordenada w[1]", "Coordenada w[2]")
  result
}

w <- c(1,1)
result <- gradientDescent(w,0.05,30000000000000000000,E,DerivadaE)
result

#resultAntiguo2b


#Ejercicio 3

# Ejercicio 1.3a

F <- function(x,y){
  (x-2)^2 + 2*(y + 2)^2 + 2*sin(2*pi*x)*sin(2*pi*y)
}

DerivadaF <- function(w){
  x = w[1]
  y = w[2]
  
  dx = 2*(2*pi*cos(2*pi*x)*sin(2*pi*y)+x-2)
  dy = 4*(pi*sin(2*pi*x)*cos(2*pi*y)+y+2)
  
  c(dx,dy)
}

gradientDescent3a <- function(w,n,maxiters,funcion,derivada){
  
  i=1
  iteraciones=0
  valores=0
  
  while(i <= maxiters ){
    
    valores[i] <- funcion(w[1],w[2])
    w = w - n * derivada(w)

    i=i+1
  }
  
  cbind(Iteraciones=1:maxiters, Valores=valores)
}

w <- c(1,1)
result <- gradientDescent3a(w,0.01,50,F,DerivadaF)
plot(result,type="b")
result <- gradientDescent3a(w,0.1,50,F,DerivadaF)
plot(result,type="b")
result

# Ejercicio 1.3b

gradientDescent3b <- function(w,n,maxiters,funcion,derivada){
  
  i=1
  coordenadas = 0
  min=3000
  
  while(i <= maxiters ){
    
    w = w - n * derivada(w)
    
    aux <- funcion(w[1],w[2])
    
    if(aux < min){
      min <- aux
      coordenadas[1] <- w[1]
      coordenadas[2] <- w[2]
    }
    
    i=i+1
  }
  
  cbind(Mínimo=min, X=coordenadas[1],Y=coordenadas[2])
}


w<-c(2.1,-2.1)
result1 <- gradientDescent3b(w,0.01,300,F,DerivadaF)
w<-c(3,-3)
result2 <- gradientDescent3b(w,0.01,300,F,DerivadaF)
w<-c(1.5,1.5)
result3 <- gradientDescent3b(w,0.01,300,F,DerivadaF)
w<-c(1,-1)
result4 <- gradientDescent3b(w,0.01,300,F,DerivadaF)

result <- rbind(result1,result2,result3,result4)
result
#---------------------

# Clasificación

#Ejercicio 2.1

# Cargar los datos:

digit.train <- read.table("../DigitosZip/zip.train", quote="\"", comment.char="", stringsAsFactors=FALSE)
digitos15.train = digit.train[digit.train$V1==1 | digit.train$V1==5,]
digitos = digitos15.train[,1]    # vector de etiquetas del train
ndigitos = nrow(digitos15.train)

grises = array(unlist(subset(digitos15.train,select=-V1)),c(ndigitos,16,16))

rm(digit.train)
rm(digitos15.train)

# Cargar Intensidad y Simetría
fsimetria <- function(A){
  A = abs(A-A[,ncol(A):1])
  -sum(A)
}

intensidadtr <- apply(grises,1,mean)

simetriatr <- apply(grises,1,fsimetria)

#Guardar datos y etiquetas (X,Y)
datosTr = cbind( apply(grises,1,mean), apply(grises,1,fsimetria))
etiquetasTr = digitos
etiquetasTr[etiquetasTr==5]=-1
etiquetasTr <- as.matrix(etiquetasTr)


# Pseudo inversa
Pseudo_inversa <- function(X, Y){
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  X <- cbind(X, 1) # cbind(X,1)
  
  x = t(X) %*% X
  pseudo = svd(x)
  aux = pseudo$v%*%diag(1/pseudo$d)%*%t(pseudo$v)
  pseudoinversa = aux%*%t(X)
  w = pseudoinversa%*%Y
  w
}

pasoARecta= function(w){
  if(length(w)!= 3)
    stop("Solo tiene sentido con 3 pesos")
  a = -w[1]/w[2]
  b = -w[3]/w[2]
  c(a,b)
}

w <- Pseudo_inversa(datosTr,etiquetasTr)

resultado_recta <- pasoARecta(w)

recta <- function(x,a,b){
  y = a*(x)+b
}

x1=0 
y1 = recta(x1,resultado_recta[1],resultado_recta[2])

x2=-5000
y2= recta(x2,resultado_recta[1],resultado_recta[2])

plot( intensidadtr, simetriatr, col=etiquetasTr+2)
#puntos_recta
lines(c(x1,x2),c(y1,y2), col="blue",type="l")


# Gradiente Descendente Estocástico MAL
SGD <- function(X, Y, eta = 0.01, umbral = 0.01, max_i = 500, seed = NULL){
  
  X = cbind(X,1)
  theta = sample(0,ncol(X),replace=TRUE)
  fin = TRUE
  
  i = 0
  stochasticList = c(1:nrow(X))
  set.seed(seed)
  while(i < max_i & fin){
    temporaryTheta = theta
    stochasticList = sample(stochasticList)
    for(i in stochasticList){
      theta = theta - eta*(-Y[i]*X[i,])/as.vector((1+exp(Y[i]*(theta%*%X[i,]))))
    }
    if(sqrt(sum((temporaryTheta - theta)^2)) < umbral)
      fin = FALSE
    i = i + 1
  }
  set.seed(NULL)
  theta
}

#######Intento 3


SGD2 <- function(X, Y, eta = 0.01, umbral = 0.1, max_i = 1000, seed = NULL, per = 0.25) {
  i = 1
  X <- cbind(X,1)
  set.seed(seed)
  dataTrain <- cbind(X,Y)
  dataTrain = dataTrain[sample(nrow(dataTrain)),]
  w = rep (0, ncol(X))
  
  while (i < max_i) {
    for (j in 1:ncol(X)) {
      w[j] = w[j] - eta * sum_SGD2 (dataTrain[1:(per * nrow(dataTrain)),], w, j)
    }
    dataTrain = dataTrain[sample(nrow(dataTrain)),]
    i = i + 1
  }
  set.seed(NULL)
  w
}

sum_SGD2 <- function (subset, w, j) {
  acc = 0.0
  indice_D = seq(1,ncol(subset)-1, by=1)
 
  for (n in 1:nrow(subset)) {
    acc = acc + subset[n,j] * (sign(sum(w*subset[n,indice_D])) - subset[n,ncol(subset)])
  }
  acc
}

w <- SGD2(datosTr,etiquetasTr,seed=1)

#Clasifica (Yprima) le da una etiqueta según la función h(w) = sign(wt*x)

Clasificar <- function(X,w){
  X <- cbind(X,1) # cbind(X,1)
  sign(X%*%w)
}

Yprima <- Clasificar(datosTr,w)
length(etiquetasTr[etiquetasTr==Yprima])

binary_error <- function(Yprima, Y){
  error <- sample(0,length(Y),replace = T)
  error[ Yprima != Y ] = 1
  error
}

Err <- function(Yprima, Y, error){
  sum(error(Yprima,Y))/length(Y)
}

Ein <- Err(Yprima,etiquetasTr,binary_error)

resultado_recta <- pasoARecta(w)

recta <- function(x,a,b){
  y = a*(x)+b
}

x1=0 
y1 = recta(x1,resultado_recta[1],resultado_recta[2])

x2=-5000
y2= recta(x2,resultado_recta[1],resultado_recta[2])

#yprima = datosTr %*% as.matrix(w[, 1:ncol(w) - 1])
#plot( simetria, intensidad, col=etiquetasTr+2, ylim=c(-250,250),xlim=c(-250,250))
plot( intensidadtr, simetriatr, col=etiquetasTr+2)
#puntos_recta
lines(c(x1,x2),c(y1,y2), col="blue",type="l")

#Ahora tendría que clasificar con la misma w las X del test
#Y calcular el Eout con la Y del test de verdad

# Cargar los datos test:

digit.test <- read.table("../DigitosZip/zip.test", quote="\"", 
                         comment.char="", stringsAsFactors=FALSE)
digitos15.test = digit.test[digit.test$V1==1 | digit.test$V1==5,]
digitos = digitos15.test[,1]    # vector de etiquetas del train
ndigitos = nrow(digitos15.test)

grises = array(unlist(subset(digitos15.test,select=-V1)),c(ndigitos,16,16))
grises = as.numeric(grises)
dim(grises)=c(49,16,16)
rm(digit.test)
rm(digitos15.test)

intensidadtst <- apply(grises,1,mean)
simetriatst <- apply(grises,1,fsimetria)

datosTest = as.matrix(cbind(intensidadtst,simetriatst))
etiquetasTest = digitos
etiquetasTest[etiquetasTest==5]=-1
etiquetasTest <- as.matrix(etiquetasTest)

Yprima <- Clasificar(datosTest,w)

Eout <- Err(Yprima, etiquetasTest,binary_error)

recta <- function(x,a,b){
  y = a*(x)+b
}

x1=0
y1 = recta(x1,resultado_recta[1],resultado_recta[2])

x2=-5000
y2= recta(x2,resultado_recta[1],resultado_recta[2])

#yprima = datosTr %*% as.matrix(w[, 1:ncol(w) - 1])
#plot( simetria, intensidad, col=etiquetasTr+2, ylim=c(-250,250),xlim=c(-250,250))
plot( intensidadtst, simetriatst, col=etiquetasTest+2)
#puntos_recta
lines(c(x1,x2),c(y1,y2), col="blue",type="l")

#########################
#2.2a
simula_unif = function (N=2,dims=2, rango = c(0,1)){
  m = matrix(runif(N*dims, min=rango[1], max=rango[2]),
             nrow = N, ncol=dims, byrow=T)
  m
}
set.seed(1)
X <- simula_unif(1000,2,c(-1,1))

plot(X)

#2.2b
f <- function(x1,x2){
  sign( (x1-0.2)^2 + x2^2 - 0.6 )
}

Y <- matrix(f(X[,1],X[,2]))

noise <- function(label, p){
  result <- label*sample(c(1, -1), size=length(label), replace=TRUE, prob= c(1-p,p))
  result
}

Y <- noise(Y,0.1)

plot(X, col=Y+2)

#2.2c

w <- SGD2(X,Y)

Yprima <- Clasificar(X,w)

Ein <- Err(Yprima,Y,binary_error)

Ein
#2.2d
vectorEi=0
vectorEo=0

for(i in 1:1000){
  X <- simula_unif(1000,2,c(-1,1))
  Y <- matrix(f(X[,1],X[,2]))
  Y <- noise(Y,0.1)
  
  w <- SGD2(X,Y,max_i = 100)
  Yprima <- Clasificar(X,w)
  vectorEi[i] <- Err(Yprima,Y,binary_error)
  
  X <- simula_unif(1000,2,c(-1,1))
  Y <- matrix(f(X[,1],X[,2]))
  Yprima <- Clasificar(X,w)
  vectorEo[i] <- Err(Yprima,Y,binary_error)
  
  print (i)
}

ValorMedio_Ein = sum(vectorEi)/length(vectorEi)
ValorMedio_Ein
ValorMedio_Eout = sum(vectorEo)/length(vectorEo)
ValorMedio_Eout
