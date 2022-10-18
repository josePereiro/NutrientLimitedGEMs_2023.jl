# # TODO: This is necessary condition but not sufficient in the context of the elementary modes of met01.
# # Document this and think about the name
# ## ------------------------------------------------------------------
# function __fill_products_connectome!(conn, S, M, N, lb, ub, met0i)
    
#     @inbounds for rxni in 1:N

#         conn[rxni] && continue # find all not-included 
            
#         if (S[met0i, rxni] < 0 && ub[rxni] > 0) # open foward consumer reactions
            
#             # include in connectome
#             conn[rxni] = true
#             # launch for each product
#             for meti in 1:M
#                 if (S[meti, rxni] > 0)
#                     # @info "found prod" rxni meti
#                     __fill_products_connectome!(conn, S, M, N, lb, ub, meti)
#                 end
#             end

#         elseif (S[met0i, rxni] > 0 && lb[rxni] < 0) # open backward consumer reactions

#             # include in connectome
#             conn[rxni] = true
#             # launch for each reactant
#             for meti in 1:M
#                 if (S[meti, rxni] < 0)
#                     __fill_products_connectome!(conn, S, M, N, lb, ub, meti)
#                 end
#             end
            
#         end

#     end # for rxni in 1:N
# end

# # Add check for bounds
# function _fill_products_connectome!(conn, S, lb, ub, met0i)
#     M, N = size(S)
#     @assert length(lb) == N && length(ub) == N && length(conn) == N
#     @assert 1 <= met0i <= M
#     return __fill_products_connectome!(conn, S, M, N, lb, ub, met0i)
# end

# """

# Fill the product connectome from `met`.
# That is, all the OPEN reactions which are and either CONSUMING `met` or any
# other metabolite that can be PRODUCED from it.
# Ex: In a network `-> A -> B <-> C`, starting from `met="A"`, it will include the middle and last reaction (Both `B` and `C` can be produced from `A` and those reactions are connecting them).

# """
# function products_connectome!(conn::Vector{Bool}, model::CoreModel, met)
#     meti = metindex(model, met)
#     _fill_products_connectome!(conn, model.S, model.xl, model.xu, meti)
#     return conn
# end
# products_connectome(model::CoreModel, met) = 
#     products_connectome!(fill(false, size(model.S, 2)), model, met)


# ## ------------------------------------------------------------------
# function __fill_precursors_connectome!(conn, S, M, N, lb, ub, met0i)
    
#     @inbounds for rxni in 1:N

#         conn[rxni] && continue # find all not-included 
            
#         if (S[met0i, rxni] > 0 && ub[rxni] > 0) # open forward producer reactions
            
#             # include in connectome
#             conn[rxni] = true
#             # launch for each "reactant"
#             for meti in 1:M
#                 if (S[meti, rxni] < 0)
#                     # @info "found prod" rxni meti
#                     __fill_precursors_connectome!(conn, S, M, N, lb, ub, meti)
#                 end
#             end

#         elseif (S[met0i, rxni] < 0 && lb[rxni] < 0) # open backward producer reactions

#             # include in connectome
#             conn[rxni] = true
#             # launch for each "product"
#             for meti in 1:M
#                 if (S[meti, rxni] > 0)
#                     __fill_precursors_connectome!(conn, S, M, N, lb, ub, meti)
#                 end
#             end
            
#         end

#     end # for rxni in 1:N
# end

# # Add check for bounds
# function _fill_precursors_connectome!(conn, S, lb, ub, met0i)
#     M, N = size(S)
#     @assert length(lb) == N && length(ub) == N && length(conn) == N
#     @assert 1 <= met0i <= M
#     return __fill_precursors_connectome!(conn, S, M, N, lb, ub, met0i)
# end

# """

# Fill the precursor connectome from `met`.
# That is, all the OPEN reactions which are either PRODUCING `met` or any
# other metabolite that is a PRECURSOR of it.
# Ex: In a network `-> A -> B <-> C`, starting from `met="C"`, it will include the middle and last reaction (Both `B` and `C` are precursors of `A` and those reactions are connecting them).

# """
# function precursors_connectome!(conn::Vector{Bool}, model::CoreModel, met)
#     meti = metindex(model, met)
#     _fill_precursors_connectome!(conn, model.S, model.xl, model.xu, meti)
#     return conn
# end
# precursors_connectome(model::CoreModel, met) = 
#     precursors_connectome!(fill(false, size(model.S, 2)), model, met)

# ## ------------------------------------------------------------------
# # All the reaction connecting from and into `met`
# function connectome(model::CoreModel, met)
#     meti = metindex(model, met)
#     conn = products_connectome(model, meti)
#     conn .|= precursors_connectome(model, meti)
#     return conn
# end

# # All the reaction connecting `met0` to `met1`
# function connectome(model::CoreModel, met0, met1)
#     conn = products_connectome(model, met0)
#     conn .&= precursors_connectome(model, met1)
#     return conn
# end