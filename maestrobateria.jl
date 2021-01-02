#Ejemplo:Problema con almacenamiento multiperiodo deterministico
#Minimiza  costo de operacion
#FUNCIONAA, BATERIA SE CARGA Y DESCARGA EN DISTINTOS PERIODOS

using JuMP
using GLPK
using Plots

# Esto hace que se puedan usar las funciones del archivo "esclavo.jl"
include("esclavo.jl")

#Demanda por periodo
demanda = [80 20 80 110] # Demanda conocida
#Definir periodos
num_periodos=4 #perido totales
num_periodos_alpha=num_periodos-1 #perido totales
periodos = 1:num_periodos

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
E_inicial=50

#Eficiencia
eta_carga=1
eta_descarga=1

#Costo por carga y descarga
Cb_carga=2
Cb_descarga=1
#Costo por perdida de carga
C_ll=200

#escenario
I_escenario=[1 1]

#DEFINIR EL MODELO ETAPA 1---------------------------------------------------
model = Model(GLPK.Optimizer)

#Define  variables de desicion
@variable(model, Pb) #Potencia bateria
@variable(model, Pg[generadores]>=0) #Potencia generadores
@variable(model, Pb_c>=0) #Potencia carga
@variable(model, Pb_d>=0) #Potencia descarga
@variable(model, Eb_final>=0) #Energia de la bateria
@variable(model, LL>=0) #Perdida de carga
@variable(model, alpha>=0) #costo futuro

#Restriccion de balance para cada etapa
@constraint(model, sum(Pg[g] for g = generadores) + Pb +LL >= demanda[1])

#Restriccion mínima y máxima capacidad generadores
for g in generadores
    @constraint(model, Pg[g] >= Pg_min[g]*I_escenario[g])
    @constraint(model, Pg[g] <= Pg_max[g]*I_escenario[g])
end

#Restriccion mínima y máxima capacidad bateria
@constraint(model, Pb >= -Pb_max)
@constraint(model, Pb <= Pb_max)

@constraint(model, Eb_final >= E_min)
@constraint(model, Eb_final <= E_max)

#No se puede cargar y descargar al mismo tiempo
@constraint(model, Pb == Pb_d-Pb_c)

#Restricciones balance de energia bateria
@constraint(model, Eb_final == E_inicial - Pb_d/eta_descarga+Pb_c*eta_carga)
@constraint(model, E_inicial >=  Pb_d/eta_descarga - Pb_c*eta_carga)

#Funcion objetivo minimo costo de operacion
@objective(model, Min, sum(Cg_var[g]*Pg[g] for g in generadores)+Pb_d*Cb_descarga+Pb_c*Cb_carga+ LL*C_ll+ alpha)


num_iteraciones = 8
iteraciones = 1:num_iteraciones

#Guarda cotas
cota_inferior_itera=zeros(1,num_iteraciones)
Cota_superior_soldef=zeros(1,num_iteraciones)

#Guardar valores periodos
#coef_pos_acumulado = zeros(1,num_periodos)
alpha_solotros=zeros(num_periodos,num_iteraciones)
Pg_solotros=zeros(num_generadores,num_iteraciones, num_periodos)
Pb_solotros=zeros(num_periodos,num_iteraciones)
Pb_csolotros=zeros(num_periodos,num_iteraciones)
Pb_dsolotros=zeros(num_periodos,num_iteraciones)
Eb_finalsolotros=zeros(num_periodos,num_iteraciones)
LL_solotros=zeros(num_periodos,num_iteraciones)

#Guardar pendiente y coef_pos corte de benders
pendiente_corte=zeros(num_periodos,num_iteraciones)
coef_pos_corte=zeros(num_periodos,num_iteraciones)

#guarda costo presente de las etapas para calcular cota superior
costo_presente_p1=0
costo_presente_etapas=zeros(num_periodos,num_iteraciones)


for i = iteraciones
    p=1
    optimize!(model)
    # La cota inferior es la solución del maestro:
    cota_inferior = objective_value(model)
    cota_inferior_itera[i]= cota_inferior
    #guardar resultado generación por tecnología p1
    Pg_solotros[:,i,p]=value.(Pg)
    Pb_solotros[p,i]=value.(Pb)
    Pb_csolotros[p,i]=value.(Pb_c)
    Pb_dsolotros[p,i]=value.(Pb_d)
    #Estado del la batería al final del p1
    Eb_finalsolp1=value.(Eb_final)
    Eb_finalsolotros[p,i]=Eb_finalsolp1
    #Perdida de carga al final del periodo 1
    LL_solotros[p,i]=value.(LL)
    #Costo futuro perido 1
    alpha_solotros[p,i]=value.(alpha)

    #Costo presente periodo 1
    costo_presente_etapas[p,i]=sum(Cg_var[g]*Pg_solotros[g,i,p] for g in generadores)+Pb_dsolotros[p,i]*Cb_descarga+Pb_csolotros[p,i]*Cb_carga+ LL_solotros[p,i]*C_ll

    for p in 2:num_periodos
        #Resuelve el subproblema periodos posteriores---------------------------
        coef_pos, pendiente, Pg_sol,  Eb_sol,LL_sol, Pb_csol, Pb_dsol,Pb_sol,alpha_sol=esclavo(p, I_escenario, Eb_finalsolotros[p-1,i], demanda[p], Eb_finalsolotros[p,:], pendiente_corte[p,:], coef_pos_corte[p,:])

        #Guardar valor corte de benders pendiente se obtiene del dual, se realiza de la siguiente manera:
        #alpha[p-1]>=coef_pos_corte[p]+pendiente_corte[p]*(Eb-Eb_sol[p-1])
        #Pendiente y coef_pos del corte a utilizar en la etapa p-1
        pendiente_corte[p-1,i]=pendiente
        coef_pos_corte[p-1,i]=coef_pos

        #Guarda costo presente de cada etapa
        costo_presente_etapas[p,i]=sum(Cg_var[g]*Pg_sol[g] for g in generadores)+Pb_dsol*Cb_descarga+Pb_csol*Cb_carga+ LL_sol*C_ll

        #guarda costo futuro del periodo
        alpha_solotros[p,i]=alpha_sol

        #guardar valores periodo
        Pg_solotros[:,i,p]=Pg_sol
        Pb_solotros[p,i]=Pb_sol
        Pb_csolotros[p,i]=Pb_csol
        Pb_dsolotros[p,i]=Pb_dsol
        Eb_finalsolotros[p,i]=Eb_sol
        LL_solotros[p,i]=LL_sol

    end
    p=1
    #corte de Benders
    @constraint(model, alpha >= coef_pos_corte[1,i] + pendiente_corte[1,i]*(Eb_final - Eb_finalsolotros[1,i]))

    # La cota superior corresponde a la suma de los costros presentes dada la solución de los esclavos:
    cota_superior = sum(costo_presente_etapas[:,i])
    Cota_superior_soldef[i]= cota_superior
    println("Cota superior: ", cota_superior)
    println("cota_inferior: ", cota_inferior)

    # Si la diferencia entre ambas cotas es muy chica, terminar el algoritmo:
    if abs(cota_superior - cota_inferior) <= 0.001
        global converge=i
        break
    end
end

#MOSTRAR SOLUCIONES
println(" ")
println("Num iteraciones: ", converge)
println("Costo Total: ", Cota_superior_soldef[converge])
println(" ")

for p in periodos
    println("Despacho PERIODO ",p,": ")
    println("Costo presente: ", costo_presente_etapas[p,converge])
    println("Gs: ", Pg_solotros[:,converge,p])
    println("Carga Bateria: ", Pb_csolotros[p,converge])
    println("Descarga Bateria: ", Pb_dsolotros[p,converge])
    println("Carga no suministrada: ", LL_solotros[p,converge])
    println(" ")
end

#Gráfico convergencia
#plot(iteraciones[1:converge],Cota_superior_soldef[1:converge], title = "Criterio de convergencia", label = "Cota superior", lw = 2)
#plot!(iteraciones[1:converge],cota_inferior_itera[1:converge],  label = "Cota inferior", lw = 2)
#xlabel!("Iteración")
#ylabel!("Valor")
#savefig("convdeter.jpg")
