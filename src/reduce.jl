function average_over_k1(k1, kbins, bk)
    Nk = length(kbins)
    bk_reduced = zeros(Nk)
    nk_reduced = zeros(Nk)
    for i in eachindex(bk)
        for j in 1:Nk
            if k1[i] >= kbins[j]
                bk_reduced[j] += bk[i]
                nk_reduced[j] += 1
            end
        end
    end
    for i in eachindex(kbins)
        bk_reduced[i] /= nk_reduced[i]
    end
    return bk_reduced
end

function average_equilateral(k1, k2, k3, max_ratio, kbins, bk)
    Nk = length(kbins)
    bk_reduced = zeros(Nk)
    nk_reduced = zeros(Nk)
    for i in eachindex(bk)
        d12 = abs((k1[i] - k2[i])/(k1[i] + k2[i]))
        d13 = abs((k1[i] - k3[i])/(k1[i] + k3[i]))
        if d12 < max_ratio && d13 < max_ratio
            for j in 1:Nk
                if k1[i] >= kbins[j]
                    bk_reduced[j] += bk[i]
                    nk_reduced[j] += 1
                end
            end
        end
    end
    for i in eachindex(kbins)
        bk_reduced[i] /= nk_reduced[i]
    end
    return bk_reduced
end

function average_squeezed(k1, k2, k3, max_ratio, kbins, bk)
    Nk = length(kbins)
    bk_reduced = zeros(Nk)
    nk_reduced = zeros(Nk)
    for i in eachindex(bk)
        d12 = abs((k1[i] - k2[i])/(k1[i] + k2[i]))
        d13 = abs((k1[i] - k3[i])/(k1[i] + k3[i]))
        if d12 < max_ratio
            for j in 1:Nk
                if k1[i] >= kbins[j]
                    bk_reduced[j] += bk[i]
                    nk_reduced[j] += 1
                end
            end
        end
    end
    for i in eachindex(kbins)
        bk_reduced[i] /= nk_reduced[i]
    end
    return bk_reduced
end

function average_folded(k1, k2, k3, max_ratio, kbins, bk)
    Nk = length(kbins)
    bk_reduced = zeros(Nk)
    nk_reduced = zeros(Nk)
    for i in eachindex(bk)
        d23 = abs((k3[i] - k2[i])/(k3[i] + k2[i]))
        if d23 < max_ratio
            for j in 1:Nk
                if k1[i] >= kbins[j]
                    bk_reduced[j] += bk[i]
                    nk_reduced[j] += 1
                end
            end
        end
    end
    for i in eachindex(kbins)
        bk_reduced[i] /= nk_reduced[i]
    end
    return bk_reduced
end