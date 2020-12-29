#Problema con almacenamiento multiperiodo estocastico
#Minimiza  costo de operacion
#resolución SDDP,masiteraciones
using JuMP
using GLPK
using CSV
using DataFrames
using Random

# Utilizar funciones del archivo "UnaEtapa.jl"
include("UnaEtapa.jl")

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
#prob=[0.8;0.1;0.1]

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

#Muestreo aleatorio
items = [1, 2, 3]
estado=zeros(Int64,num_periodos)
#for j in 1:40
#    for i = 1:length(estado)
#        estado[i] = sample(items, Weights(prob))
#    end
#    println(estado)
#end
Esc_Resueltos=[1 1 1;1 1 2;1 1 3;1 2 1;1 2 2;1 2 3;1 3 1; 1 3 2;1 3 3]
I=zeros(num_generadores,num_periodos)
Resultados_pg=zeros(9,num_periodos,num_generadores)
Resultados_pb=zeros(9,num_periodos)
Resultados_eb=zeros(9,num_periodos)
Resultados_ll=zeros(9,num_periodos)

#e=9
for e in 1:9
        rama= Esc_Resueltos[e,:]
        for p in 1:num_periodos
                I[:,p]=I_realizaciones[rama[p],:]
        end
        eb_final_sol, pg_sol, pb_sol, ll_sol=UnaEtapa(I, demanda, num_periodos,e_in)
        for p in 1:num_periodos
                Resultados_pg[e,p,:]=pg_sol[:,p]
                Resultados_pb[e,p]=pb_sol[p]
                Resultados_eb[e,p]=eb_final_sol[p]
                Resultados_ll[e,p]=ll_sol[p]
        end
end
