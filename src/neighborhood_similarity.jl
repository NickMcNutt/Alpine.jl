using VantagePointTrees

import Base.LinAlg.BLAS: gemm!
import VantagePointTrees: initneighbors!, updateneighbors!, knn!, sumsqrtgt, Metric

function initneighbors!(clusters::Clusters, index::Int)
    @inbounds begin
        nd = clusters.dist
        k = clusters.num_neighbors
        @simd for j in 1:k
            nd[j, index] = Inf
        end
    end
end

function updateneighbors!{T<:Real}(clusters::Clusters, index::Int, d::T, i::Int)
    @inbounds begin
        nd = clusters.dist
        ni = clusters.indices_neighbors
        k = clusters.num_neighbors

        j = 1
        while d > nd[j, index]
            j += 1
        end

        for m in k:-1:j+1
            nd[m, index] = nd[m - 1, index]
            ni[m, index] = ni[m - 1, index]
        end

        nd[j, index] = d
        ni[j, index] = i

        τ = nd[k, index]
        return τ
    end
end

function knn!{T<:Real}(clusters::Clusters, index::Int, τ::T, metric::PeriodicSqEuclidean{T}, points, node::Node, x::T, y::T, z::T)
    @inbounds begin
        i = node.index
        i == 0 && return τ

        d::T = evaluate(metric, points, i, x, y, z)::T

        if d < τ && d > 0.001
            τ = updateneighbors!(clusters, index, d, i)
        end

        node.isleaf && return τ

        μ = node.distance
        # if τ + d <= μ
        if !sumsqrtgt(τ, d, μ)
            τ = knn!(clusters, index, τ, metric, points, node.left, x, y, z)
        else
            τ = knn!(clusters, index, τ, metric, points, node.right, x, y, z)
            # if d <= τ + μ
            if sumsqrtgt(τ, μ, d)
                τ = knn!(clusters, index, τ, metric, points, node.left, x, y, z)
            end
        end

        return τ
    end
end

function knn!{T<:Real}(clusters::Clusters, index::Int, τ::T, metric::PeriodicSqEuclidean{T}, points, node::Node, p::Int)
    @inbounds begin
        i = node.index
        i == 0 && return τ

        d::T = evaluate(metric, points, i, p)::T

        if d < τ && d > 0.001
            τ = updateneighbors!(clusters, index, d, i)
        end

        node.isleaf && return τ

        μ = node.distance
        # if τ + d <= μ
        if !sumsqrtgt(τ, d, μ)
            τ = knn!(clusters, index, τ, metric, points, node.left, p)
        else
            τ = knn!(clusters, index, τ, metric, points, node.right, p)
            # if d <= τ + μ
            if sumsqrtgt(τ, μ, d)
                τ = knn!(clusters, index, τ, metric, points, node.left, p)
            end
        end

        return τ
    end
end

function knn!{T<:Real}(clusters::Clusters, index::Int, tree::VPTree, x::T, y::T, z::T)
    initneighbors!(clusters, index)
    knn!(clusters, index, Inf, tree.metric, tree.points, tree.root, x, y, z)
end

function knn!(clusters::Clusters, index::Int, tree::VPTree, p::Int)
    initneighbors!(clusters, index)
    @inbounds clusters.indices_central[index] = p
    knn!(clusters, index, Inf, tree.metric, tree.points, tree.root, p)
    @inbounds distances!(tree.metric.bw, clusters.coords, index, tree.points, clusters.indices_neighbors, p)
end

