
#####################
# Plotting functions
#####################

myColors = [[0.9,0.9,0.1],
            [0.9,0.1,0.9],
            [0.1,0.9,0.1],
            [0.1,0.1,0.9]]

def plotSimpleMisinc(df, letterorder=['C', 'A', 'T', 'G'], colors=myColors, ):
    misinc_df, lb_df, ub_df = simpleMisincStats(df)
    x = misinc_df.index