module JoubertSyndrome

using DataFrames
using CSV
using DataStructures: SortedDict, Reverse
using StatsBase
using CategoricalArrays
using HypothesisTests
using Plots
using StatsPlots
using FreqTables


DIR_OUT = "/Users/olivierlabayle/Documents/JoubertSyndrome"

DATAPATH = joinpath(DIR_OUT, "AA_these_joubert_descriptif_862021_VF.csv")


PHEN_OF_INTEREST = ("Neonatal breathing dysregulations",
                    "Epilepsy",
                    "Obesity or overweight",
                    "Macrocephaly",
                    "Retinal dystrophy",
                    "Coloboma",
                    "Kidney defect",
                    "Liver defect",
                    "Orofacial anomalies",
                    "Multiple pulmonary infections",
                    "Endocrinal and/or genital features ",
                    "polydactyly",
                    "pure JS")

GENES_OF_INTEREST = ("CC2D2A",
                    "CPLANE1",
                    "AHI1",
                    "CEP290",
                    "TMEM67",
                    "KIAA0586",
                    "INPP5E",
                    "CSPP1",
                    "RPGRIP1L",
                    "OFD1",
                    "MKS1",
                    "NPHP1")


data = DataFrame(CSV.File(DATAPATH; missingstrings=["na"]))


# ###########################################################################
#         PVALUES                                                           #
# ###########################################################################

function getfigures(x, y)
    ct = freqtable(x, y)
    odds_phen_if_gen = (ct[2,2]/ct[2,1]) / (ct[1,2]/ct[1,1])
    se = sqrt(sum(1 ./ ct))
    lower_bound = exp(log(odds_phen_if_gen) - 1.96*se)
    upper_bound = exp(log(odds_phen_if_gen) + 1.96*se)
    min_freq = minimum(ct)
    reschi = ChisqTest(ct)
    resfish = FisherExactTest(ct[1,1], ct[1,2], ct[2,1], ct[2,2])
    return min_freq, reschi, resfish, odds_phen_if_gen, lower_bound, upper_bound
end


function gene_phenotype_it(data::DataFrame, gene, phenotype)
    data_restrict = dropmissing(data[:, ["GENE", phenotype]])
    x = categorical(data_restrict.GENE .== gene)
    y = categorical(data_restrict[:, phenotype] .== "oui" )
    return getfigures(x, y)
end


function phenotype_phenotype_it(data::DataFrame, phenotype1, phenotype2)
    data_restrict = dropmissing(data[:, [phenotype1, phenotype2]])
    x = categorical(data_restrict[:, phenotype1] .== "oui")
    y = categorical(data_restrict[:, phenotype2] .== "oui")
    return getfigures(x, y)
end


function gene_phenotype_it_results(data::DataFrame, genes, phenotypes)
    results = DataFrame(Gene = String[], Phenotype = String[], 
                PVAL_CHI2 = Float64[], PVAL_FISH = Float64[], MIN_FREQ = Int[],
                OR = Float64[], OR_LB = Float64[], OR_UB = Float64[], 
                RES_CHI2 = Any[], RES_FISH = Any[])
    for gene in genes
        for phenotype in phenotypes
            min_freq, reschi, resfish, odds_phen_if_gen, lower_bound, upper_bound = gene_phenotype_it(data, gene, phenotype)
            push!(results, 
            (gene, phenotype, pvalue(reschi), pvalue(resfish), min_freq, odds_phen_if_gen, lower_bound, upper_bound, reschi, resfish))
        end
    end
    return results
end

function phenotype_phenotype_it_results(data::DataFrame, phenotypes)
    results = DataFrame(Phenotype1 = String[], Phenotype2 = String[], 
                        PVAL_CHI2 = Float64[], PVAL_FISH = Float64[], 
                        MIN_FREQ = Int[], OR = Float64[], OR_LB = Float64[], OR_UB = Float64[],
                        RES_CHI2 = Any[], RES_FISH = Any[])
    for phenotype1 in phenotypes
        for phenotype2 in phenotypes
            if phenotype1 != phenotype2
                min_freq, reschi, resfish, odds_phen_if_gen, lower_bound, upper_bound = phenotype_phenotype_it(data, phenotype1, phenotype2)
                push!(results, 
                (phenotype1, phenotype2, pvalue(reschi), pvalue(resfish), min_freq, odds_phen_if_gen, lower_bound, upper_bound, reschi, resfish))
            end
        end
    end
    return results
end


gene_phen_results = gene_phenotype_it_results(data, GENES_OF_INTEREST, PHEN_OF_INTEREST)

phen_phen_results = phenotype_phenotype_it_results(data, PHEN_OF_INTEREST)

CSV.write(joinpath(DIR_OUT, "gene_phen.csv"), gene_phen_results[:, Not(["RES_CHI2", "RES_FISH"])])
CSV.write(joinpath(DIR_OUT, "phen_phen.csv"), phen_phen_results[:, Not(["RES_CHI2", "RES_FISH"])])


# ###########################################################################
#         PLOTS                                                             #
# ###########################################################################


genes_to_plot = ["CC2D2A", "CPLANE1", "AHI1", "CEP290", "TMEM67", "OFD1", "NPHP1", "CSPP1"]

phens_to_plot = ["Neonatal breathing dysregulations", "Epilepsy", "Retinal dystrophy", "Coloboma",
"Kidney defect", "Liver defect", "Orofacial anomalies", "polydactyly", "Multiple pulmonary infections",
"Endocrinal and/or genital features ", "pure JS"]


cols_to_xticks_fr = Dict(
    "Multiple pulmonary infections" => "Infections respiratoires \nrépétées",
    "polydactyly" => "Polydactylie",
    "Endocrinal and/or genital features " => "Anomalies \nendocrino-génitales",
    "Orofacial anomalies" => "Anomalies orofaciales",
    "Liver defect" => "Atteinte hépatique",
    "Kidney defect" =>  "Atteinte rénale",
    "Retinal dystrophy" => "Atteinte rétinienne",
    "Coloboma" => "Colobome",
    "Epilepsy" => "Epilepsie",
    "Neonatal breathing dysregulations" => "Troubles respiratoires \nnéonataux",
    "pure JS" => """JS \"pur\""""
 )
 

 cols_to_xticks_en = Dict(
    "Multiple pulmonary infections" => "Multiple pulmonary \ninfections",
    "polydactyly" => "Polydactyly",
    "Endocrinal and/or genital features " => "Endocrinal and/or \ngenital features",
    "Orofacial anomalies" => "Orofacial anomalies",
    "Liver defect" => "Liver disease",
    "Kidney defect" =>  "Kidney disease",
    "Retinal dystrophy" => "Retinal disease",
    "Coloboma" => "Coloboma",
    "Epilepsy" => "Epilepsy",
    "Neonatal breathing dysregulations" => "Neonatal breathing \ndysregulations",
    "pure JS" => """\"Pure\" JS"""
    )


function myfreq(x)
    nomiss_x = skipmissing(x)
    nbyes = sum(nomiss_x .== "oui")
    nb = length(collect(nomiss_x))
    return round(100*nbyes/nb, digits=2)
end

function getstars(gene_phen_results::DataFrame, gene, phen)
    pvals = filter(row->(row.Gene == gene && row.Phenotype == phen),  
            gene_phen_results)[1, [:PVAL_FISH, :PVAL_CHI2]]
    pvals = Vector(pvals)
    if all(pvals .< 0.001)
        return " (***)"
    elseif all(pvals .< 0.01)
        return " (**)"
    elseif all(pvals .< 0.05)
        return " (*)"
    else
        return ""
    end
end

function Base.unique(ctg::CategoricalArray)
    l = levels(ctg)
    newctg = CategoricalArray(l)
    levels!(newctg, l)
end

function save_plots(data, genes_to_plot, phens_to_plot, gene_phen_results, cols_to_xticks;
                     dirout=DIR_OUT,
                     lang="en",
                     titlebase="Phenotype affection for gene ", 
                     xlabel="Phenotypes", 
                     ylabel="Percentage of affected individuals",
                     groups=["Mutated", "Not mutated"])
    export_df = DataFrame()
    for gene in genes_to_plot
        data[!, gene] = data[!, "GENE"] .== gene
        grouped = groupby(data, gene)
        freqs = combine(grouped, phens_to_plot .=> myfreq)
        mutation_true_freqs = Vector(filter(row->row[gene] === true, freqs)[1, Not(gene)])
        mutation_false_freqs = Vector(filter(row->row[gene] === false, freqs)[1, Not(gene)])
        xticks = [cols_to_xticks[phen] * getstars(gene_phen_results, gene, phen) for phen in phens_to_plot]
        xticks = categorical(repeat(xticks, outer=2), levels=xticks)
        p = groupedbar(xticks, [mutation_true_freqs mutation_false_freqs], 
                    bar_position=:dodge, 
                    bar_width=0.7,
                    # xlabel=xlabel,
                    xtickfontsize=7,
                    title=titlebase * gene,
                    rotation=45,
                    fontfamily="Times",
                    ylims=(0, 100),
                    legend=:best,
                    ylabel=ylabel,
                    group=repeat(groups, inner=length(phens_to_plot)))
        savefig(p, joinpath(dirout, "$(gene)_barplot_$lang.png"))
        export_df[!, "$(gene)_mutation_true_freqs"] = mutation_true_freqs
        export_df[!, "$(gene)_mutation_false_freqs"] = mutation_false_freqs
    end
    return export_df
end

data_to_plot = data[:, vcat(phens_to_plot, ["GENE"])]

save_plots(data_to_plot, genes_to_plot, phens_to_plot, gene_phen_results, cols_to_xticks_en;
                     dirout=DIR_OUT,
                     lang="en",
                     titlebase="Phenotype affection for gene ", 
                     xlabel="Phenotypes", 
                     ylabel="Percentage of affected individuals",
                     groups=["Mutated", "Not mutated"])

export_df = save_plots(data_to_plot, genes_to_plot, phens_to_plot, gene_phen_results, cols_to_xticks_fr;
                     dirout=DIR_OUT,
                     lang="fr",
                     titlebase="Phénotype du gène ", 
                     xlabel="Phénotypes", 
                     ylabel="Pourcentage d'individus affectés",
                     groups=["Muté", "Non muté"])

# export_df[!, "Phenotypes"] = phens_to_plot
# CSV.write(joinpath(DIR_OUT, "freqs.csv"), export_df)


end
