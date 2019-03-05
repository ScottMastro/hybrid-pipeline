import plot_nucmer_helper as plotter
import matplotlib.pyplot as plt

name="slc9a3_insert"
startPos = 205903000

#prefix="C:/Users/Scott/Desktop/plotnucmer/export"
prefix="~/export"
files=[prefix + "/CF066B2D.contigs.fasta.PILON2.fasta.slc9a3_insert.coords"
       #prefix + "/nfixed_362.fasta"  + ".coords",
      # prefix + "/CF062B2D.canu.fasta." + name + ".coords",
      # prefix + "/OSK7121.supernova2.fasta." + name + ".coords"
        # "~/export/nfixed.fa." + name + ".coords"
       ]

snp_files=[    
        prefix + "/CF066B2D.contigs.fasta.PILON2.fasta.slc9a3_insert.snps"
        #prefix + "/nfixed_362.fasta" + ".snps",
         #  prefix + "/CF062B2D.canu.fasta." + name + ".snps",
      #prefix + "/OSK7121.supernova2.fasta." + name + ".snps",
       #"~/export/nfixed.fa." + name + ".snps"
      ]

tracknames = ["Hybrid Assembly", "PacBio Assembly (mixed haplotypes)", "10x Assembly (one haplotype)"]
gff_file=prefix + "/" + name + ".gff"
source="ensembl"

lenFilter = 1000
idFilter = 0.92

def plot_all():

    fig, ax = plt.subplots(1)
    totalwidth = plotter.get_totalwidth(files)
    plotter.initialize_plot(fig, ax, tracknames, totalwidth)
    plotter.add_genes(ax, gff_file, totalwidth)

    for i, file in enumerate(files):
        print("reading alignments for " + file)
        alignments=plotter.read_alignments(file, lenFilter=lenFilter, idFilter=idFilter)
        plotter.add_alignments(ax, alignments, i, len(files), totalwidth,
                               snp_file=snp_files[i])
                    
    fig.canvas.draw()
    plotter.adjust_xaxis(ax, startPos)
        
    plt.savefig("region.svg")
    plt.show()
    
plot_all()
