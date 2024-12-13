= Methods

== Segmentatino

== Alignment to events
We aligned the trace of each segmented cell to the occurrence of poke events. 

== Decoder
To construct the decoder, we first took the average activity in a 1 second window leading up to each poke event. We then subtracted the average activity of this activity across all poke events, for each cell. This resulted a $n_("cells") x n_("events")$ activity matrix $X$. To reduce the number of dimensions that went into the decoder, we performed singular value decomposition (SVD) on this matrix. That is, we found $U$, $Sigma$ and $V$ such that

$X = U Sigma V^t$,

where $U$ and $V$ are matrices with orthonormal columns and $Sigma$ is a diagonal matrix containing the singular values of $X$. We formed a lower dimensional basis for the subspace by retaining enough singular vectors to account for 75% of the singular values, and concatenating these vectors into a projection matrix $W$. We then projected the activity matrix $X$ onto the subspace spanned by $W$, creating a lower dimensional activity matrix $Y = W X$. 
To quantify how much information about the poke direction (left vs right) was contained in the subspace, we fit a linear discriminant function to distinguish between activity preceeding a left poke from activity preceeding a right poke, resulting in another projection matrix $L$. To get a time-resolved readout of the poke direction information, we went back to the raw cell traces. We subtracted the mean across time for each trace, and then projected the matrix containing the traces for all cells across time onto the subspace spanned by $W$, and then onto the linear discriminant function represented by the matrix $L$, to give a one dimensional readout $y$,

$y = L W X$,

where a postiive value represented evidence for a left poke and negative values evidence for a right poke. The further away from zero the value of $y$ is at any point in time, the more confident the decoder was about its choice.

