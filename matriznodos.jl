#Problema con almacenamiento multiperiodo estocastico
#Minimiza  costo de operacion
#Resolución directa, formulación por nodos
#Definitiva
using JuMP
using GLPK
using CSV
using DataFrames

#PARAMETROS CARACTERIZACION DEL ARBOL------------------------------------------
#definir nodos
nodos_finales=9
numero_niveles=3
orden=3
#nodos_totales=(numero_niveles*nodos_finales-1)/(numero_niveles-1)
nodos_totales=13
nodos=1:nodos_totales
#definir periodes
periodos=zeros(Int64,nodos_totales)
for n in nodos
    if n==1
        periodos[n]=1
    elseif n==2 ||n==3 ||n==4
        periodos[n]=2
    else
        periodos[n]=3
    end
end

#Definir vectores de ceros para almacenar padres y hermanos
padre_matrix=zeros(Int64,nodos_totales,nodos_totales) # 1 si el nodo x es padre del nodo y
padre_numero=zeros(Int64,nodos_totales) #numero de nodo padre
hermanos=zeros(Int64,nodos_totales,nodos_totales) # 1 si el nodo x es hernado del nodo y

#DATOS DEMANDA----------------------------------------------------------
demanda = [80 20 80] # Demanda conocida
demanda_nodo =zeros(nodos_totales)

for n in nodos
    demanda_nodo[n]=demanda[periodos[n]]
end
println(demanda_nodo)
#DATOS CENTRALES TERMICAS------------------------------------------------
#Numero de centrales generadores
num_generadores=2
generadores=1:num_generadores
#Costo operacion
cg_vr=[40 100]

#Limites tecnicos
pg_max=[50 50]
pg_min=[0 0]

#DATOS BATERIA-----------------------------------------------------------------
#Limites tecnicos
pb_max=50
pb_min=0
e_max=50
e_min=0
e_in=50 #estado de carga inicial

#Eficiencia
eta_ch=0.9 #carga
eta_di=1 #descarga

#Costo por carga y descarga
cb_ch=5 #carga
cb_di=5 #descarga

#PERDIDA DE CARGA-----------------------------------------------
VoLL=1200


#PROBABILIDADES DE LOS NODOS---------------------------------------------------
prob=zeros(nodos_totales) #vector para almacenar la provabilidad de los nodos
I_nodos=zeros(Int64,nodos_totales,2)
tasa_falla_an=5 #occ/años
tasa_falla_h=5/8760 #occ/hora
prob_ramas=[1-tasa_falla_h*2,tasa_falla_h,tasa_falla_h]
I_ramas=[1 1; 1 0; 0 1]

#definir matriz de padres
for p in nodos, h in nodos
    if h ==2*p+(p-1) || h ==2*p+(p-1)+1 || h ==2*p+(p-1)+2
        padre_matrix[p,h]=1
        padre_numero[h]=p
    end
end

#Definir matriz de hermanos
for x in nodos, y in nodos
    if padre_numero[x] == padre_numero[y]
        hermanos[x,y]=1
    end
end

#calculo probabilidades de cada nodo
prob[1]=1
for n in 2:orden:nodos_totales
    prob[n]=prob[padre_numero[n]]*prob_ramas[1]
    prob[n+1]=prob[padre_numero[n+1]]*prob_ramas[2]
    prob[n+2]=prob[padre_numero[n+2]]*prob_ramas[3]
end
println(prob)

#calculo disponibilidad de cada nodo
prob[1]=1
for n in 2:orden:nodos_totales
    prob[n]=prob[padre_numero[n]]*prob_ramas[1]
    prob[n+1]=prob[padre_numero[n+1]]*prob_ramas[2]
    prob[n+2]=prob[padre_numero[n+2]]*prob_ramas[3]
end
#calculo disponibilidad en cada nodo
I_nodos[1,:]=I_ramas[1,:]

for n in 2:orden:nodos_totales
    I_nodos[n,:]=I_ramas[1,:]
    I_nodos[n+1,:]=I_ramas[2,:]
    I_nodos[n+2,:]=I_ramas[3,:]
end
println(I_nodos)

#DEFINIR EL MODELO ETAPA 1---------------------------------------------------
model = Model(GLPK.Optimizer)

#Define  variables de desicion estapa 1
@variable(model, pb[nodos]) #Potencia bateria
@variable(model, pg[nodos,generadores]>=0) #Potencia generadores
@variable(model, pb_c[nodos]>=0) #Potencia carga
@variable(model, pb_d[nodos]>=0) #Potencia descarga
@variable(model, eb_final[nodos]>=0) #Energia de la bateria
@variable(model, ll[nodos]>=0) #Perdida de carga

#Restricciones etapa 1------------------------------------------------------
#Restriccion de balance
for n in nodos
    @constraint(model, sum(pg[n,g] for g = generadores) + pb[n] +ll[n] >= demanda_nodo[n])
end
#Restriccion mínima y máxima capacidad generadores
for n in nodos, g in generadores
    @constraint(model, pg[n,g] >= pg_min[g]*I_nodos[n,g])
    @constraint(model, pg[n,g] <= pg_max[g]*I_nodos[n,g])
end

#Restriccion mínima y máxima capacidad bateria
for n in nodos
    @constraint(model, pb[n] >= -pb_max)
    @constraint(model, pb[n] <= pb_max)
end

for n in nodos
    @constraint(model, eb_final[n] >= e_min)
    @constraint(model, eb_final[n] <= e_max)
end

#No se puede cargar y descargar al mismo tiempo
for n in nodos
    @constraint(model, pb[n] == pb_d[n]-pb_c[n])
end

#Restricciones balance de energia bateria nodo 1
@constraint(model, eb_final[1] == e_in - pb_d[1]/eta_di +pb_c[1]*eta_ch)
@constraint(model, e_in >=  pb_d[1]/eta_di - pb_c[1]*eta_ch)

#Restricciones balance de energia bateria nodo 2 en adelante
for n in 2:nodos_totales
    @constraint(model, eb_final[n] == eb_final[padre_numero[n]] - pb_d[n]/eta_di +pb_c[n]*eta_ch)
    @constraint(model, eb_final[padre_numero[n]] >=  pb_d[n]/eta_di - pb_c[n]*eta_ch)
end

#Funcion objetivo minimo costo de operacion
@objective(model, Min, sum(prob[n]*(sum(cg_vr[g]*pg[n,g] for g in generadores)+pb_d[n]*cb_di+pb_c[n]*cb_ch+ ll[n]*VoLL) for n in nodos))
optimize!(model)

#SOLUCIONES POR NODO
costo = objective_value(model)
pg_sol=value.(pg)
pb_sol=value.(pb)
pb_c_sol=value.(pb_c)
pb_d_sol=value.(pb_d)
eb_sol=value.(eb_final)
ll_sol=value.(ll)


#MOSTRAR SOLUCIONES
println("Costo total: ", costo)
println(" ")

resultados = DataFrame(Gs1 = pg_sol[:,1],
               Gs2 = pg_sol[:,2],
               Carga_bateria = pb_c_sol,
               Descarga_bateria = pb_d_sol,
               Estado_carga_finp = eb_sol,
               Perdida_carga = ll_sol,
               )
CSV.write("Resultadosnodos.csv", resultados)
