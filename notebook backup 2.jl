### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ e9c95740-c126-11ee-2377-8d36b4e39b09
begin
	using Plots, Images
end

# ╔═╡ b1e9ef52-de4e-4023-8805-4bacf8152f3b
begin 
	# First, create a random NxN matrix 
	function initial_state(N::Int64)
		return rand(N+2,N+2) .< 0.5
	end

	# plot matrix as image
	function img_mtx(M)
		return Gray.(M)
	end
end

# ╔═╡ dc1ab299-9187-45be-82ae-cafbf5c79b31
begin
	N = 100
	M = initial_state(N)
	img_mtx(M)
end

# ╔═╡ 1d740e01-b1c0-4cbf-9a4d-30f1e1566d2b
Int8(9)

# ╔═╡ b3a85ec5-261f-4c34-ae9b-9aac730264da
begin 
	# Calculate the total energy 
	# of a given distribution of 2D spins 
	function cspin(x::Bool)::Int8
		return Int8(x)*2-1
	end
	function sumcspin(i::Int,j::Int)::Int
		return cspin(M[i,j+1])+cspin(M[i,j-1])+cspin(M[i+1,j])+cspin(M[i-1,j])
	end
	function E_ij(J::Float64,M,i::Int,j::Int)::Float64
		# Given i,j 
		# get the sum of their
		# Neigbors contributions to the total
		# Energy, like this:
		#    	 (i,j+1)
		# (i-1,j)(i ,j )(i+1,j)
		# 		(i ,j-1)
		S = cspin(!M[i,j])*sumcspin(i,j)
		return -J*S
	end
end

# ╔═╡ fff49d55-9625-4956-86a8-330359d1f41e
begin 
	function E_tot(M,J)
		# Get Matrix size
		N = size(M,1)-2 # remember to subtract two
		S = 0.0 
		for i in 2:N
			for j in 2:N
				S += E_ij(J,M,i,j)
			end
		end
		return S
	end
end

# ╔═╡ da215597-d8fc-4d55-aaa6-f99d89a716fb
E_old = E_tot(M,0.1)

# ╔═╡ 1f749e66-d8e9-4bb8-b392-819559326efc
begin
	# Now we will create a function that will 
	function update_monte_carlo(M,J,β,E_old)
		# Ensure periodicity
		#M[begin,:] = M[end,:]
		#M[:,begin] = M[:,end]
		# Get Matrix size
		N = size(M,1)-2 # remember to subtract two
		# Select a random site
		#i,j = rand(2:N,2)
		E_new = E_old
		for i in 2:N
			for j in 2:N
				# Get its nearest neighbors
				# Energy contribution to its change
				E_new = E_old + E_ij(J,M,i,j)
				Pratio = exp(β*(E_old-E_new))
				if Pratio ≥ rand() 
					M[i,j] = !M[i,j]
					E_old = E_new
				end
			end 
		end
		return [M,E_new] # Aceitar a mudança
	end
end

# ╔═╡ 9f8725eb-1877-40be-8d6e-2862d3e636aa
begin
	function update_(M,J::Float64,β::Float64,n::Int64)
		E_old = E_tot(M,J)
		for i in 1:n
			M,E_old = update_monte_carlo(M,J,β,E_old)
			#println(E_old)
		end
		return M
	end
end

# ╔═╡ b1d6ad08-4f2d-4cf6-a165-e22b51bfb5d9
img_mtx(M)

# ╔═╡ 7415e815-cd55-4d93-b497-829824de671a
img_mtx(update_(M,0.03,100.0,1))

# ╔═╡ 9286d595-af8c-4975-993a-6758c2e07216
md"# Need to fix periodic boundary conditions"

# ╔═╡ Cell order:
# ╠═e9c95740-c126-11ee-2377-8d36b4e39b09
# ╠═b1e9ef52-de4e-4023-8805-4bacf8152f3b
# ╠═dc1ab299-9187-45be-82ae-cafbf5c79b31
# ╠═1d740e01-b1c0-4cbf-9a4d-30f1e1566d2b
# ╠═b3a85ec5-261f-4c34-ae9b-9aac730264da
# ╠═fff49d55-9625-4956-86a8-330359d1f41e
# ╠═da215597-d8fc-4d55-aaa6-f99d89a716fb
# ╠═1f749e66-d8e9-4bb8-b392-819559326efc
# ╠═9f8725eb-1877-40be-8d6e-2862d3e636aa
# ╠═b1d6ad08-4f2d-4cf6-a165-e22b51bfb5d9
# ╠═7415e815-cd55-4d93-b497-829824de671a
# ╠═9286d595-af8c-4975-993a-6758c2e07216
