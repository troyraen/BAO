# save wtheta_list to a file

# convert to pandas df
def wtheta_list_to_df(wtheta_list):
    """Assumes wtheta_list has nm (# mocks) elements of form [bcens, wtheta(bcen)].
    Returns dataframe with nb (# bins) rows, nm+2 columns.
        Column 1 is bcens; 2 thru nm+1 are wtheta(bcen) for each mock,
        nm+2 is average wtheta of all mocks."""
    nb = len(bcens)
    nm = len(wtheta_list)

    wtheta_df = pd.DataFrame()
    bcens0 = wtheta_list[0][0]
    for i, mock in enumerate(wtheta_list):
        np.testing.assert_array_equal(bcens0, mock[0], \ # ensure bcens[i] == bcens0
            err_msg='bcens0 != bcens{}\n(also check floating point rounding error)'.format(i))
        wtheta_df['mock{}_wtheta'.format(i)] = pd.Series(mock[1]) # add column of wtheta

    wtheta_df["wtheta_avg"] = wtheta_df.mean(axis=1) # average the mocks
    wtheta_df.insert(0, "bin_center", bcens0) # add bin_centers to front

    return wtheta_df


# wtheta to file write to file
wtheta_df = wtheta_list_to_df(wtheta_list)
__, colnames = wtheta_df.axes()
fname = './wtheta.data'
np.savetxt(fname, wtheta_df.to_numpy(), delimiter=',', header=','.join(colnames))
