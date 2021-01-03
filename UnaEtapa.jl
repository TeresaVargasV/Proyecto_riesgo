#Problema con almacenamiento multiperiodo estocastico
#Minimiza  costo de operacion
#Problema de una etapa
using JuMP
using GLPK
using CSV
using DataFrames


#DATOS CENTRALES TERMICAS------------------------------------------------
#Numero de centrales generadores
#num_generadores=2
num_generadores=5
generadores=1:num_generadores
#Costo operacion
#cg_vr=[40 100]
cg_vr=[5 12 25 15 20]


#Limites tecnicos
#pg_max=[50 50]
#pg_min=[0 0]
pg_max=[400 350 300 250 200]
pg_min=[0 0 0 0 0]

#DATOS BATERIA-----------------------------------------------------------------
#Limites tecnicos
#pb_max=50
#pb_min=0
#e_max=50
#e_min=0
pb_max=150
pb_min=0
e_max=150
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
function UnaEtapa(I_escenario, demanda, num_periodos,eb_in)
    #DEFINIR EL MODELO ETAPA 1---------------------------------------------------
    model = Model(GLPK.Optimizer)

    #Define  variables de desicion estapa
    @variable(model, pb[1:num_periodos]) #Potencia bateria
    @variable(model, pg[generadores,1:num_periodos]>=0) #Potencia generadores
    @variable(model, pb_c[1:num_periodos]>=0) #Potencia carga
    @variable(model, pb_d[1:num_periodos]>=0) #Potencia descarga
    @variable(model, eb_final[1:num_periodos]>=0) #Energia final de la bateria
    #@variable(model, eb_inicial[num_periodos]>=0) #Energia inicial de la bateria
    @variable(model, ll[1:num_periodos]>=0) #Perdida de carga
    #@variable(model, alpha[num_periodos]>=0) #funcion de costo futuro

    #Restricciones etapa------------------------------------------------------
    
    for t in 1:num_periodos
        #Restriccion de balance
        @constraint(model, sum(pg[g,t] for g = generadores) + pb[t] +ll[t] >= demanda[t])

        #Restriccion mínima y máxima capacidad generadores
        for g in generadores
            @constraint(model, pg[g,t] >= pg_min[g]*I_escenario[g,t])
            @constraint(model, pg[g,t] <= pg_max[g]*I_escenario[g,t])
        end

        #Restriccion mínima y máxima capacidad bateria
        @constraint(model, pb[t] >= -pb_max)
        @constraint(model, pb[t] <= pb_max)

        #retriccion capacidad embalse
        @constraint(model, eb_final[t] >= e_min)
        @constraint(model, eb_final[t] <= e_max)

        #entrega potencia de salida
        @constraint(model, pb[t] == pb_d[t]-pb_c[t])
    end

    #Restricciones de inventario periodo 1
    @constraint(model, eb_final[1] == eb_in - pb_d[1]/eta_di +pb_c[1]*eta_ch)

    #Restricciones de inventario a partir de periodo 2
    for t in 2:num_periodos
        @constraint(model, eb_final[t] == eb_final[t-1] - pb_d[t]/eta_di +pb_c[t]*eta_ch)
    end

    #Funcion objetivo minimo costo de operacion
    @objective(model, Min, sum(sum(cg_vr[g]*pg[g,t] for g in generadores)+pb_d[t]*cb_di+pb_c[t]*cb_ch+ ll[t]*VoLL for t in 1:num_periodos))
    optimize!(model)

    #Obtener soluciones para los cortes
    #coef_pos = objective_value(model)
    #pendiente = dual.(fijar)
    eb_final_sol=value.(eb_final)
    pg_sol=value.(pg)
    pb_sol=value.(pb)
    #alpha_sol=value.(alpha)
    ll_sol=value.(ll)

    # devolver lo necesario para los cortes
    return eb_final_sol, pg_sol, pb_sol, ll_sol
end
