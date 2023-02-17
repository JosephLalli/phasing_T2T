#!/usr/bin/env python

def plot_ideogram(UCSC_cytoband_coords, datapoints, 
                  data_height=0.4, data_padding=0.1,
                  figsize=(6, 8),
                  chromosome_list = ['chr%s' % i for i in list(range(1, 23)) + ['M', 'X', 'Y']],
                  inclusive=False):
    """
    Slightly modified from https://gist.github.com/daler/c98fc410282d7570efc3.
    Demonstrates plotting chromosome ideograms and genes (or any features, really)
    using matplotlib.

    1) Assumes a file from UCSC's Table Browser from the "cytoBandIdeo" table,
    saved as "ideogram.txt". Lines look like this::

        #chrom  chromStart  chromEnd  name    gieStain
        chr1    0           2300000   p36.33  gneg
        chr1    2300000     5300000   p36.32  gpos25
        chr1    5300000     7100000   p36.31  gneg

    2) Assumes another dataframe containing ['chrom', 'start', 'end', 'name']
       columns. These will be plotted along the chromosomes.    

    """

    from matplotlib import pyplot as plt
    from matplotlib.collections import BrokenBarHCollection
    import pandas as pd


    # Here's the function that we'll call for each dataframe (once for chromosome
    # ideograms, once for genes).  The rest of this script will be prepping data
    # for input to this function
    #
    def chromosome_collections(df, y_positions, height,  **kwargs):
        """

        Yields BrokenBarHCollection of features that can be added to an Axes
        object.

        Parameters
        ----------

        df : pandas.DataFrame
            Must at least have columns ['chrom', 'start', 'end', 'color']. If no
            column 'width', it will be calculated from start/end.

        y_positions : dict
            Keys are chromosomes, values are y-value at which to anchor the
            BrokenBarHCollection

        height : float
            Height of each BrokenBarHCollection

        Additional kwargs are passed to BrokenBarHCollection
        """
        del_width = False
        if 'width' not in df.columns:
            del_width = True
            df['width'] = df['end'] - df['start']
        for chrom, group in df.groupby('chrom'):
            yrange = (y_positions[chrom], height)
            xranges = group[['start', 'width']].values
            yield BrokenBarHCollection(
                xranges, yrange, facecolors=group['colors'], **kwargs)
        if del_width:
            del df['width']


    # Height of each ideogram
    chrom_height = 1

    # Height of the gene track. Should be smaller than `chrom_spacing` in order to
    # fit correctly
    data_height = data_height

    # Spacing between consecutive ideograms
    chrom_spacing = data_height * 2 + data_padding * 2


    # Keep track of the y positions for ideograms and genes for each chromosome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    data_ybase = {}
    chrom_centers = {}

    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        data_ybase[chrom] = ybase - data_height - data_padding
        ybase += chrom_height + chrom_spacing

    # Read in ideogram.txt, downloaded from UCSC Table Browser
    ideo = pd.read_table(
        UCSC_cytoband_coords,
        skiprows=1,
        names=['chrom', 'start', 'end', 'name', 'gieStain']
    )

    # Filter out chromosomes not in our list
    ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]

    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    # Colors for different chromosome stains
    color_lookup = {
        'gneg': (1., 1., 1.),
        'gpos25': (.6, .6, .6),
        'gpos50': (.4, .4, .4),
        'gpos75': (.2, .2, .2),
        'gpos100': (0., 0., 0.),
        'acen': (.8, .4, .4),
        'gvar': (.8, .8, .8),
        'stalk': (.9, .9, .9),
    }

    # Add a new column for colors
    ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])


    # Same thing for genes
    datapoints = datapoints[datapoints.chrom.apply(lambda x: x in chromosome_list)]
    # if whole chromosome to be displayed, add start and end data points
    if inclusive:
        data_list=list()
        for chrom, df in datapoints.groupby('chrom'):
            # print(df)
            empty_line=pd.DataFrame([{'chrom':chrom, 'start':-1, 'end':-1, 'name':-1}])
            chrom_end = ideo.loc[ideo.chrom==chrom, 'end'].max()
            if df.iat[0, 1] != 0:
                df = pd.concat((empty_line, df))
                df.iat[0, 1] = 0
                df.iat[0, 2] = df.iat[1, 1]
            if df.iat[0, 2] != chrom_end:
                df = pd.concat((df, empty_line))
                df.iat[-1, 1] = df.iat[-2, 2]
                df.iat[-1, 2] = chrom_end
            data_list.append(df)
        datapoints = pd.concat(data_list)


    datapoints['width'] = datapoints.end - datapoints.start
    datapoints['colors'] = '#2243a8'
    datapoints.loc[datapoints.index %2 == 0, 'colors'] = '#8c5939'


    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Now all we have to do is call our function for the ideogram data...
    print("adding ideograms...")
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
        ax.add_collection(collection)

    # ...and the gene data
    print("adding genes...")
    for collection in chromosome_collections(
        datapoints, data_ybase, data_height, alpha=0.5, linewidths=0
    ):
        ax.add_collection(collection)

    # Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.axis('tight')
    return ax
