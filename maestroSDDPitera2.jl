#Problema con almacenamiento multiperiodo estocastico
#Minimiza  costo de operacion
#resolución SDDP,masiteraciones
using JuMP
using GLPK
using CSV
using DataFrames
using Random

# Utilizar funciones del archivo "esclavoSDDP.jl"
include("esclavoSDDPm2.jl")

#DEFINIR PERIODOS Y ESTADOS----------------------------------------
num_periodos=3 #perido totales
periodos = 1:num_periodos
num_real=3
realizaciones=1:num_real

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
I_realizaciones=[1 1;1 0; 0 1]
prob=[1-5/8760*2;5/8760;5/8760]
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
@variable(model, pb) #Potencia bateria
@variable(model, pg[generadores]>=0) #Potencia generadores
@variable(model, pb_c>=0) #Potencia carga
@variable(model, pb_d>=0) #Potencia descarga
@variable(model, eb_final>=0) #Energia de la bateria
@variable(model, ll>=0) #Perdida de carga
@variable(model, alpha>=0) #funcion de costo futuro

#Restricciones etapa 1------------------------------------------------------
#Restriccion de balance
@constraint(model, sum(pg[g] for g = generadores) + pb +ll >= demanda[1])

#Restriccion mínima y máxima capacidad generadores
for g in generadores
    @constraint(model, pg[g] >= pg_min[g])
    @constraint(model, pg[g] <= pg_max[g])
end

#Restriccion mínima y máxima capacidad bateria
@constraint(model, pb >= -pb_max)
@constraint(model, pb <= pb_max)

@constraint(model, eb_final >= e_min)
@constraint(model, eb_final <= e_max)

#No se puede cargar y descargar al mismo tiempo
@constraint(model, pb == pb_d-pb_c)

#Restricciones balance de energia bateria
@constraint(model, eb_final == e_in - pb_d/eta_di +pb_c*eta_ch)
@constraint(model, e_in >=  pb_d/eta_di - pb_c*eta_ch)

#Funcion objetivo minimo costo de operacion
@objective(model, Min, sum(cg_vr[g]*pg[g] for g in generadores)+pb_d*cb_di+pb_c*cb_ch+ ll*VoLL+ alpha)

#DEFINIR NUMERO DE ITERACIONES-----------------------------------------------
num_iteraciones=15
iteraciones = 1:num_iteraciones
# Muestreo
num_muestreo=2
muestreo=zeros(Int64,num_muestreo)

#Guardar pendiente y coef_pos corte de benders
eb_final_p=zeros(num_periodos,num_muestreo,num_iteraciones)
pendiente_cortes=zeros(num_periodos,num_muestreo,num_iteraciones)
coef_cortes=zeros(num_periodos,num_muestreo,num_iteraciones)
pendiente_corte_aux=zeros(num_real,num_muestreo,num_periodos)
coef_corte_aux=zeros(num_real,num_muestreo,num_periodos)
#Guardar costo presente, y calculo de cotas para de cada iteracion
costo_presente_aux=zeros(num_periodos,num_muestreo,num_iteraciones)
z_max=zeros(num_iteraciones)
z_min=zeros(num_iteraciones)
var_z=0
converge=0
ms1=0
ms2=0
for i = iteraciones
    optimize!(model)
    p=1
    #Obtener soluciones periodo 1
    costo_maestro = objective_value(model)
    for m in 1:num_muestreo
        eb_final_p[p,m,i]=value.(eb_final)
    end
    alpha_m_sol=value.(alpha)
    z_min[i]=costo_maestro
    for m in 1:num_muestreo
        costo_presente_aux[p,m,i]=costo_maestro-alpha_m_sol
    end

    #Etapa forward
    for p in 2:num_periodos
        global ms1=1
        global ms2=1
        while ms1==ms2
            global ms1=rand((1,2,3))
            global ms2=rand((1,2,3))
        end
        muestreo[1]=ms1
        muestreo[2]=ms2
        #Resuelve el subproblema periodos posteriores muestreo 1---------------------------
        coef_pos, pendiente, eb_final_sol, alpha_sol,pg_per_sol,pb_per_sol=esclavoSDDPm2(I_realizaciones[muestreo[1],:], eb_final_p[p-1,1,i], demanda[p],coef_cortes[p,:,:], pendiente_cortes[p,:,:], eb_final_p[p,:,:],num_iteraciones,num_muestreo)

        #Guardar valor de la variable con acople temporal
        eb_final_p[p,1,i]=eb_final_sol
        costo_presente_aux[p,1,i]=coef_pos-alpha_sol

        #Resuelve el subproblema periodos posteriores muestreo 2---------------------------
        coef_pos, pendiente, eb_final_sol, alpha_sol,pg_per_sol,pb_per_sol=esclavoSDDPm2(I_realizaciones[muestreo[2],:], eb_final_p[p-1,2,i], demanda[p],coef_cortes[p,:,:], pendiente_cortes[p,:,:], eb_final_p[p,:,:],num_iteraciones,num_muestreo)

        #Guardar valor de la variable con acople temporal
        eb_final_p[p,2,i]=eb_final_sol
        costo_presente_aux[p,2,i]=coef_pos-alpha_sol
    end

    #Ver convergencia
    z_max[i]=costo_presente_aux[1,1,i]+ (sum(costo_presente_aux[p,e,i] for p in 2:num_periodos, e in 1:num_muestreo))*1/num_muestreo
    #global var_z=sqrt(1/(num_muestreo^2)*sum((z_max[i]-sum(costo_presente_aux[p,m,i] for p in 1:num_periodos))^2 for m in 1:num_muestreo))
    if z_max[i] - z_min[i]<=10
        global converge=i
        break
    end

    #backward
    for p in num_periodos:-1:2
        for m in 1:num_muestreo
            for e in realizaciones
                #Resuelve el subproblema periodos posteriores---------------------------
                coef_pos, pendiente, eb_final_sol, alpha_sol,pg_per_sol,pb_per_sol=esclavoSDDPm2(I_realizaciones[e,:,:], eb_final_p[p-1,m,i], demanda[p],coef_cortes[p,:,:], pendiente_cortes[p,:,:], eb_final_p[p,:,:],num_iteraciones,num_muestreo)
                pendiente_corte_aux[e,m,p-1]=pendiente
                coef_corte_aux[e,m,p-1]=coef_pos
            end
            pendiente_cortes[p-1,m,i]=sum(pendiente_corte_aux[e,m,p-1]*prob[e] for e in realizaciones)
            coef_cortes[p-1,m,i]=sum(coef_corte_aux[e,m,p-1]*prob[e] for e in realizaciones)
        end
    end
    p=1
    #corte de Benders periodo 1
    for m in 1:num_muestreo
        @constraint(model, alpha >= coef_cortes[1,m,i] + pendiente_cortes[1,m,i]*(eb_final - eb_final_p[1,m,i]))
    end
end

#SOLUCIONES PERIODO 1
#se resuelve de nuevo la función objetivo
optimize!(model)
costo_maestro = objective_value(model)
pg_m_sol=value.(pg)
pb_m_sol=value.(pb)
pb_c_m_sol=value.(pb_c)
pb_d_m_sol=value.(pb_d)
eb_final_m_sol=value.(eb_final)
ll_m_sol=value.(ll)
alpha_m_sol=value.(alpha)

println(" ")
println("Num iteraciones: ", converge)
println("Costo Total: ", costo_maestro)
println(" ")
for g in generadores
    println("Despacho Gs ",g," :", pg_m_sol[g])
end
println("Carga bateria: ", pb_c_m_sol)
println("Descarga bateria: ", pb_d_m_sol)
println("Estado de carga al final del periodo: ", eb_final_m_sol)
println("Perdida de carga: ", ll_m_sol)
println("Funcion de costo futuro: ", alpha_m_sol)
println(" ")
