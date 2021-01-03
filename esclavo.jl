#Problema con almacenamiento multiperiodo deterministico
#Minimiza  costo de operacion
#Ejemplo bateria, FUNCIONA

using JuMP
using GLPK


#DATOS CENTRALES TERMICAS------------------------------------------------
#Costo operacion
Cg_var=[40 100]

#Numero de centrales generadores
num_generadores=2
generadores=1:num_generadores

#Limites tecnicos
Pg_max=[50 50]
Pg_min=[0 0]

#DATOS BATERIA-----------------------------------------------------------------
#Limites tecnicos
Pb_max=50
Pb_min=0
E_max=50
E_min=0

#Eficiencia
eta_carga=1
eta_descarga=1

#Costo por carga y descarga
Cb_carga=2
Cb_descarga=1

#Costo por perdida de carga
C_ll=1000



function esclavo(p, I_escenario, E_anterior, demanda,Eb_result, pendiente_dual, objetivo)
    #DEFINIR EL MODELO----------------------------------------------------------
    model = Model(GLPK.Optimizer)

    #Define  variables de desicion
    @variable(model, Pb) #Potencia bateria
    @variable(model, Pg[generadores]>=0) #Potencia generadores
    @variable(model, Pb_c>=0) #Potencia carga
    @variable(model, Pb_d>=0) #Potencia descarga
    @variable(model, Eb>=0) #Energia de la bateria
    @variable(model, LL>=0) #Perdida de carga
    @variable(model, E_tmenos1>=0) #Energia al final de la etapa enterior
    @variable(model, alpha>=0) #costo futuro

    #Restriccion de balance para cada etapa
    @constraint(model, sum(Pg[g] for g = generadores) + Pb +LL == demanda)

    #Restriccion mínima y máxima capacidad generadores
    for g in generadores
        @constraint(model, Pg[g] >= Pg_min[g]*I_escenario[g])
        @constraint(model, Pg[g] <= Pg_max[g]*I_escenario[g])
    end

    #Restriccion mínima y máxima capacidad bateria
    @constraint(model, Pb >= -Pb_max)
    @constraint(model, Pb <= Pb_max)
    @constraint(model, Eb >= E_min)
    @constraint(model, Eb <= E_max)

    #Balance carga y descarga
    @constraint(model, Pb == Pb_d-Pb_c)

    #Restricciones balance de energia bateria
    @constraint(model, Eb == E_tmenos1 - Pb_d/eta_descarga+Pb_c*eta_carga)
    @constraint(model, E_tmenos1 >=  Pb_d/eta_descarga - Pb_c*eta_carga)

    #fijar energia almacenada inicial
    @constraint(model, fijar, E_tmenos1 == E_anterior)

    #Corte de Benders de la forma alpha[p-1]>=coef_pos_corte[p]+pendiente_corte[p]*(Eb-Eb_resut[p-1])
    for i in 1:length(Eb_result)
        @constraint(model, alpha >= objetivo[i] + pendiente_dual[i]*(Eb - Eb_result[i]))
    end

    #Funcion objetivo minimo costo de operacion
    @objective(model, Min, sum(Cg_var[g]*Pg[g] for g in generadores)+Pb_d*Cb_descarga+Pb_c*Cb_carga+ LL*C_ll+ alpha)

    optimize!(model)

    coef_pos = objective_value(model)
    pendiente = dual.(fijar)
    Pg_sol=value.(Pg)
    Pb_sol=value.(Pb)
    Pb_csol=value.(Pb_c)
    Pb_dsol=value.(Pb_d)
    Eb_sol=value.(Eb)
    LL_sol=value.(LL)
    alpha_sol=value.(alpha)
    return coef_pos, pendiente, Pg_sol, Eb_sol,LL_sol, Pb_csol, Pb_dsol,Pb_sol,alpha_sol
end
