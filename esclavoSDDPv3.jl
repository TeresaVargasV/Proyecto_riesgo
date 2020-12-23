#Problema con almacenamiento multiperiodoestocastico con SDDP
#Minimiza  costo de operacion
#Problema esclavo con un muestreo de 2
using JuMP
using GLPK
using CSV
using DataFrames

#DEFINIR PERIODOS Y ESTADOS----------------------------------------
num_periodos=8 #perido totales
periodos = 1:num_periodos

#DATOS DEMANDA----------------------------------------------------------
demanda = [80 20 60 50 20 150 200 180] # Demanda conocida

#DATOS CENTRALES TERMICAS------------------------------------------------
#Numero de centrales generadores
num_generadores=5
generadores=1:num_generadores
#Costo operacion
cg_vr=[40 100 60 55 120]

#Limites tecnicos
pg_max=[50 50 50 50 50]
pg_min=[0 0 0 0 0]

#DATOS BATERIA-----------------------------------------------------------------
#Limites tecnicos
pb_max=50
pb_min=0
e_max=50
e_min=0
#e_in=50 #estado de carga inicial

#Eficiencia
eta_ch=0.9 #carga
eta_di=1 #descarga

#Costo por carga y descarga
cb_ch=5 #carga
cb_di=5 #descarga

#PERDIDA DE CARGA-----------------------------------------------
VoLL=1200
function esclavoSDDPv3(I_escenario, e_anterior, demanda, coef_cortes, pendiente_cortes, eb_result,num_itera,num_muestreo)
    #DEFINIR EL MODELO ETAPA 1---------------------------------------------------
    model = Model(GLPK.Optimizer)

    #Define  variables de desicion estapa
    @variable(model, pb) #Potencia bateria
    @variable(model, pg[generadores]>=0) #Potencia generadores
    @variable(model, pb_c>=0) #Potencia carga
    @variable(model, pb_d>=0) #Potencia descarga
    @variable(model, eb_final>=0) #Energia final de la bateria
    @variable(model, eb_inicial>=0) #Energia inicial de la bateria
    @variable(model, ll>=0) #Perdida de carga
    @variable(model, alpha>=0) #funcion de costo futuro

    #Restricciones etapa------------------------------------------------------
    #Restriccion de balance
    @constraint(model, sum(pg[g] for g = generadores) + pb +ll >= demanda)

    #Restriccion mínima y máxima capacidad generadores
    for g in generadores
        @constraint(model, pg[g] >= pg_min[g]*I_escenario[g])
        @constraint(model, pg[g] <= pg_max[g]*I_escenario[g])
    end

    #Restriccion mínima y máxima capacidad bateria
    @constraint(model, pb >= -pb_max)
    @constraint(model, pb <= pb_max)

    @constraint(model, eb_final >= e_min)
    @constraint(model, eb_final <= e_max)

    #No se puede cargar y descargar al mismo tiempo
    @constraint(model, pb == pb_d-pb_c)

    #restriccion para obtener dual
    @constraint(model, fijar, eb_inicial == e_anterior)

    #Restricciones balance de energia bateria
    @constraint(model, eb_final == eb_inicial - pb_d/eta_di +pb_c*eta_ch)
    @constraint(model, eb_inicial >=  pb_d/eta_di - pb_c*eta_ch)

    #Corte de Benders de la forma alpha[p-1]>=coef_pos_corte[p]+pendiente_corte[p]*(Eb-Eb_resut[p-1])

    for  m in 1:num_muestreo, i in 1:num_itera
        @constraint(model, alpha >= coef_cortes[m,i] + pendiente_cortes[m,i]*(eb_final - eb_result[m,i]))
    end


    #Funcion objetivo minimo costo de operacion
    @objective(model, Min, sum(cg_vr[g]*pg[g] for g in generadores)+pb_d*cb_di+pb_c*cb_ch+ ll*VoLL+ alpha)
    optimize!(model)
    #println("iteracion del esclavo", i)
    #Obtener soluciones para los cortes
    coef_pos = objective_value(model)
    pendiente = dual.(fijar)
    eb_final_sol=value.(eb_final)
    pg_sol=value.(pg)
    pb_sol=value.(pb)
    alpha_sol=value.(alpha)
    println("test, ",value.(ll))

    # devolver lo necesario para los cortes
    return coef_pos, pendiente, eb_final_sol, alpha_sol,pg_sol,pb_sol
end
