#Problema con almacenamiento multiperiodo deterministico
#Minimiza  costo de operacion
#resolución directa
using JuMP
using GLPK
using CSV
using DataFrames

# Leer los archivos con datos
periodo2_file = CSV.read("Prob_t2.csv", delim=';')
periodo3_file = CSV.read("Prob_t3.csv", delim=';')

# Obtener las columnas de los archivos (nombre[!,:Nombrecolumna] )
prob_t2 = periodo2_file[!,:Prob_falla] #probabilidad de transision t2
prob_t3 = periodo3_file[!,:Prob_falla] #probabilidad de transision t3
Ig1_t2 = periodo2_file[!,:G1] #estado fallado no fallado t2, g1
Ig2_t2 = periodo2_file[!,:G2] #estado fallado no fallado t2, g2
Ig1_t3 = periodo3_file[!,:G1] #estado fallado no fallado t3, g1
Ig2_t3 = periodo3_file[!,:G2] #estado fallado no fallado t3, g2
padre_t3 = periodo3_file[!,:nodo_padre1] #nodo de origen

#DEFINIR PERIODOS Y ESTADOS----------------------------------------
num_periodos=3 #perido totales
periodos = 1:num_periodos
num_nod_p2=length(prob_t2)
nod_p2 = 1:num_nod_p2
num_nod_p3=length(prob_t3)
nod_p3 = 1:num_nod_p3

#DATOS DEMANDA----------------------------------------------------------
demanda = [80 20 80] # Demanda conocida

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

#DEFINIR EL MODELO ETAPA 1---------------------------------------------------
model = Model(GLPK.Optimizer)

#Define  variables de desicion estapa 1
@variable(model, pb_p1) #Potencia bateria
@variable(model, pg_p1[generadores]>=0) #Potencia generadores
@variable(model, pb_c_p1>=0) #Potencia carga
@variable(model, pb_d_p1>=0) #Potencia descarga
@variable(model, eb_final_p1>=0) #Energia de la bateria
@variable(model, ll_p1>=0) #Perdida de carga

#Define  variables de desicion estapa 2
@variable(model, pb_p2[nod_p2]) #Potencia bateria
@variable(model, pg_p2[generadores,nod_p2]>=0) #Potencia generadores
@variable(model, pb_c_p2[nod_p2]>=0) #Potencia carga
@variable(model, pb_d_p2[nod_p2]>=0) #Potencia descarga
@variable(model, eb_final_p2[nod_p2]>=0) #Energia de la bateria
@variable(model, ll_p2[nod_p2]>=0) #Perdida de carga

#Define  variables de desicion estapa 3
@variable(model, pb_p3[nod_p3]) #Potencia bateria
@variable(model, pg_p3[generadores,nod_p3]>=0) #Potencia generadores
@variable(model, pb_c_p3[nod_p3]>=0) #Potencia carga
@variable(model, pb_d_p3[nod_p3]>=0) #Potencia descarga
@variable(model, eb_final_p3[nod_p3]>=0) #Energia de la bateria
@variable(model, ll_p3[nod_p3]>=0) #Perdida de carga

#Restricciones etapa 1------------------------------------------------------
#Restriccion de balance
@constraint(model, sum(pg_p1[g] for g = generadores) + pb_p1 +ll_p1 >= demanda[1])

#Restriccion mínima y máxima capacidad generadores
for g in generadores
    @constraint(model, pg_p1[g] >= pg_min[g])
    @constraint(model, pg_p1[g] <= pg_max[g])
end

#Restriccion mínima y máxima capacidad bateria
@constraint(model, pb_p1 >= -pb_max)
@constraint(model, pb_p1 <= pb_max)

@constraint(model, eb_final_p1 >= e_min)
@constraint(model, eb_final_p1 <= e_max)

#No se puede cargar y descargar al mismo tiempo
@constraint(model, pb_p1 == pb_d_p1-pb_c_p1)

#Restricciones balance de energia bateria
@constraint(model, eb_final_p1 == e_in - pb_d_p1/eta_di +pb_c_p1*eta_ch)
@constraint(model, e_in >=  pb_d_p1/eta_di - pb_c_p1*eta_ch)

#Restricciones etapa 2-------------------------------------------------------------
#Restriccion de balance
for n in nod_p2
    @constraint(model, sum(pg_p2[g,n] for g = generadores) + pb_p2[n] +ll_p2[n] >= demanda[2])
end
#Restriccion mínima y máxima capacidad generadores
for n in nod_p2
    @constraint(model, pg_p2[1,n] >= pg_min[1]*Ig1_t2[n])
    @constraint(model, pg_p2[1,n] <= pg_max[1]*Ig1_t2[n])
    @constraint(model, pg_p2[2,n] >= pg_min[2]*Ig2_t2[n])
    @constraint(model, pg_p2[2,n] <= pg_max[2]*Ig2_t2[n])
end


for n in nod_p2
    #Restriccion mínima y máxima capacidad bateria
    @constraint(model, pb_p2[n] >= -pb_max)
    @constraint(model, pb_p2[n] <= pb_max)

    @constraint(model, eb_final_p2[n] >= e_min)
    @constraint(model, eb_final_p2[n] <= e_max)
    #No se puede cargar y descargar al mismo tiempo
    @constraint(model, pb_p2[n] == pb_d_p2[n]-pb_c_p2[n])
end

for n in nod_p2
    #Restricciones balance de energia bateria
    @constraint(model, eb_final_p2[n] == eb_final_p1 - pb_d_p2[n]/eta_di +pb_c_p2[n]*eta_ch)
    @constraint(model, eb_final_p1 >=  pb_d_p2[n]/eta_di - pb_c_p2[n]*eta_ch)
end

#Restricciones etapa 3-------------------------------------------------------------
#Restriccion de balance
for n in nod_p3
    @constraint(model, sum(pg_p3[g,n] for g = generadores) + pb_p3[n] +ll_p3[n] >= demanda[3])
end
#Restriccion mínima y máxima capacidad generadores
for n in nod_p3
    @constraint(model, pg_p3[1,n] >= pg_min[1]*Ig1_t3[n])
    @constraint(model, pg_p3[1,n] <= pg_max[1]*Ig1_t3[n])
    @constraint(model, pg_p3[2,n] >= pg_min[2]*Ig2_t3[n])
    @constraint(model, pg_p3[2,n] <= pg_max[2]*Ig2_t3[n])
end


for n in nod_p3
    #Restriccion mínima y máxima capacidad bateria
    @constraint(model, pb_p3[n] >= -pb_max)
    @constraint(model, pb_p3[n] <= pb_max)

    @constraint(model, eb_final_p3[n] >= e_min)
    @constraint(model, eb_final_p3[n] <= e_max)
    #No se puede cargar y descargar al mismo tiempo
    @constraint(model, pb_p3[n] == pb_d_p3[n]-pb_c_p3[n])
end

for n in nod_p3
    #Restricciones balance de energia bateria
    @constraint(model, eb_final_p3[n] == eb_final_p2[padre_t3[n]] - pb_d_p3[n]/eta_di +pb_c_p3[n]*eta_ch)
    @constraint(model, eb_final_p2[padre_t3[n]] >=  pb_d_p3[n]/eta_di - pb_c_p3[n]*eta_ch)
end


#Funcion objetivo minimo costo de operacion
@objective(model, Min, sum(cg_vr[g]*pg_p1[g] for g in generadores)+pb_d_p1*cb_di+pb_c_p1*cb_ch+ ll_p1*VoLL+sum((sum(cg_vr[g]*pg_p2[g,n] for g in generadores)+ pb_d_p2[n]*cb_di+pb_c_p2[n]*cb_ch+ ll_p2[n]*VoLL)*prob_t2[n] for n in nod_p2)+ sum((sum(cg_vr[g]*pg_p3[g,n] for g in generadores)+ pb_d_p3[n]*cb_di+pb_c_p3[n]*cb_ch+ ll_p3[n]*VoLL)*prob_t3[n] for n in nod_p3))
optimize!(model)

#Obtener soluciones periodo 1
costo_p1 = objective_value(model)
pg_p1_sol=value.(pg_p1)
pb_p1_sol=value.(pb_p1)
pb_c_p1_sol=value.(pb_c_p1)
pb_d_p1_sol=value.(pb_d_p1)
eb_p1_sol=value.(eb_final_p1)
ll_p1_sol=value.(ll_p1)

#Obtener soluciones periodo 2
pg_p2_sol=value.(pg_p2)
pb_p2_sol=value.(pb_p2)
pb_c_p2_sol=value.(pb_c_p2)
pb_d_p2_sol=value.(pb_d_p2)
eb_p2_sol=value.(eb_final_p2)
ll_p2_sol=value.(ll_p2)

#Obtener soluciones periodo 3
pg_p3_sol=value.(pg_p3)
pb_p3_sol=value.(pb_p3)
pb_c_p3_sol=value.(pb_c_p3)
pb_d_p3_sol=value.(pb_d_p3)
eb_p3_sol=value.(eb_final_p3)
ll_p3_sol=value.(ll_p3)

#MOSTRAR SOLUCIONES
println("Costo total: ", costo_p1)
println(" ")

#SOLUCIONES PERIODO 1
for g in generadores
    println("Despacho Gs ",g," :", pg_p1_sol[g])
end
println("Carga bateria: ", pb_c_p1_sol)
println("Descarga bateria: ", pb_d_p1_sol)
println("Estado de carga al final del periodo: ", eb_p1_sol)
println("Perdida de carga: ", ll_p1_sol)
println(" ")

#SOLUCIONES PERIODO 2
for n in nod_p2
    for g in generadores
        println("Despacho Gs ",g," escenario ",n," :", pg_p2_sol[g,n])
    end
    println("Carga bateria escenario ",n," :", pb_c_p2_sol[n])
    println("Descarga bateria escenario ",n," :", pb_d_p2_sol[n])
    println("Estado de carga al final del periodo 2 escenario ",n," :", eb_p2_sol[n])
    println("Perdida de carga escenario ",n," :", ll_p2_sol[n])
    println(" ")
end
#SOLUCIONES PERIODO 3
for n in nod_p3
    for g in generadores
        println("Despacho Gs ",g," escenario ",n," :", pg_p3_sol[g,n])
    end
    println("Carga bateria escenario ",n," :", pb_c_p3_sol[n])
    println("Descarga bateria escenario ",n," :", pb_d_p3_sol[n])
    println("Estado de carga al final del periodo 2 escenario ",n," :", eb_p3_sol[n])
    println("Perdida de carga escenario ",n," :", ll_p3_sol[n])
    println(" ")
end
a=["AKANKSHA", "TANYA", "PREETIKA", "VRINDA", "JAHNVI"]
pb_c_p3_sol=value.(pb_c_p3)
pb_d_p3_sol=value.(pb_d_p3)
eb_p3_sol=value.(eb_final_p3)
ll_p3_sol=value.(ll_p3)
# Creando DataFrames
resultados_p1 = DataFrame(Equipo = ["Gs1","Gs2", "Carga_bateria", "Descarga_bateria", "Estado_carga_finp","Perdida_carga"],
               Valores = [pg_p1_sol[1], pg_p1_sol[2], pb_c_p1_sol, pb_d_p1_sol,eb_p1_sol,ll_p1_sol],
               )
resultados_p2 = DataFrame(Gs1 = pg_p2_sol[1,:],
               Gs2 = pg_p2_sol[2,:],
               Carga_bateria = pb_c_p2_sol,
               Descarga_bateria = pb_d_p2_sol,
               Estado_carga_finp = eb_p2_sol,
               Perdida_carga = ll_p2_sol,
               )
resultados_p3 = DataFrame(Gs1 = pg_p3_sol[1,:],
               Gs2 = pg_p3_sol[2,:],
               Carga_bateria = pb_c_p3_sol,
               Descarga_bateria = pb_d_p3_sol,
               Estado_carga_finp = eb_p3_sol,
               Perdida_carga = ll_p3_sol,
               )
CSV.write("ResultadosP1.csv", resultados_p1)
CSV.write("ResultadosP2.csv", resultados_p2)
CSV.write("ResultadosP3.csv", resultados_p3)
