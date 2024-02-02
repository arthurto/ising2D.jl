# Definindo uma estrutura de dados que vai ser um bool 
mutable struct Site
    a::Bool
end 

function get_value(site::Site)::Int8
    return (true => +1 ,false => -1)[site.a]
end

function sum(sites::Vector{Site})::Int
    return sum(get_value.(sites))
end

abstract type Lattice <: Matrix{Site} end

function lattice(N::Int)::Lattice 
    return [Site(rand(Bool)) for _ in 1:N, _ in 1:N]
end

function index(i::Int64,nc::Int64)::Int64
    return mod(i,nr)+1
end 

function sum_neighbors(lattice::Lattice,i::Int,j::Int)::Int8 
    nr,nc = size(lattice)
    # Com certeza deve ter um jeito mais elegante de escrever isso 
    # Mas essa foi a melhor forma que eu encontrei de somar os vizinhos respeitando 
    # as condições periódicas de contorno (topologia de um torus)
    s = sum([
        lattice[index(i+1,nr)+1,index(j,nc)],
        lattice[index(i-1,nr),index(j,nc)],
        lattice[index(i,nr),index(j+1,nc)],
        lattice[index(i,nr),index(j-1,nc)]
        ])
    return s
end

function sum_neighbors(J::Real,lattice::Lattice,i::Int,j::Int)::Int8
    # Summing the neighbors of i,j
    return -J*get_value(lattice[index(i),index(j)])*sum_neighbors(lattice,i,j)
end

function sum_neighbors(J::Real,lattice::Lattice,i::Int,j::Int,flip=true)::Int8
    # Summing the neighbors of i,j, with i,j flipped
    return -J*get_value(!lattice[index(i),index(j)])*sum_neighbors(lattice,i,j)
end

function sum_neighbors(J::Real,lattice::Latttice)::Real
    # Total energy
    nr,nc = size(lattice)
    S = Real(0)
    for i in 1:nr,j in 1:nc
        S += sum_neighbors(J,lattice,i,j)
    end
    return S
end

function sum_neighbors(J::Real,lattice::Latttice,flip=true)::Real
    # Total energy for flipped
    nr,nc = size(lattice)
    S = Real(0)
    for i in 1:nr,j in 1:nc
        S += sum_neighbors(J,lattice,i,j,true)
    end
    return S
end

function flip!(lattice::Lattice,i::Int,j::Int)
    # Flipping a spin in the i,j position,
    # or, more generally, changing the current state 
    # to some other (in this case, only one) 
    # different state.
    # In this case we're using booleans
    # so the opposite of a boolean is its negation:
    # !boolean 
    lattice[index(i),index(j)] = !lattice[index(i),index(j)]
end

function evolve_state!(J::Real,β::Real,lattice::Lattice,i::Int,j::Int,E_old::Float64)
    E_new = sum_neighbors(J,lattice,i,j,true)
    
    # A variação de energia 
    # vai determinar se a alteração 
    # vai ser realizada ou não
    ΔE = E_new-E_old
    # Se a energia diminuiu, estão a troca é aceita 
    if ΔE < 0 
        flip!(lattice,i,j)
        E_old = E_new 
    else 
        r = rand() # Gera um número aleatório 
        P = exp(-β*ΔE) # Calcular o peso de Boltzmann para a transição 
        # E_old -> E_new no caso em que E_new > E_old 
        if P ≥ r # Se o peso estatístico de Boltzmann for 
            # maior que o número real aleatório r, teremos 
            # uma flip também
            flip!(lattice,i,j)
            E_old = E_new
        end
    end 
end