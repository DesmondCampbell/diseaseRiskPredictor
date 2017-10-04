# replace missings with their transpose elements (which should be non-missing)
fnCompleteSymmetricMatrix <- function( mSymm )
{
  mbx <- is.na(mSymm)
  #XXXX NB this is a GOTCHA! mSymm[mbx] <- mSymm[t(mbx)] #doesn't produce a symmetric matrix
  mSymm[mbx] <- t(mSymm)[(mbx)]

  # validate
  if(any(is.na(mSymm))) stop("Symmetric matrix contains NAs")

  mSymm
}

# creates a matrix
fnInitMatrix <- function( vsRow, vsCol, fillValue=NA)
{
  mMatrix <- matrix(fillValue,nrow=length(vsRow),ncol=length(vsCol))
  rownames(mMatrix) <- vsRow
  colnames(mMatrix) <- vsCol
  mMatrix
}

# copies source matrix into destination matrix where row/colnames match
fnCopyIntoMatrix <- function(mDest,mSrc)
{
  vbxRow <- rownames(mDest) %in% rownames(mSrc)
  vbxCol <- colnames(mDest) %in% colnames(mSrc)
  if( !sum(vbxRow)| !sum(vbxCol)) return(mDest)
  vsRow <- rownames(mDest)[vbxRow]
  vsCol <- colnames(mDest)[vbxCol]
  mDest[vsRow,vsCol] <- mSrc[vsRow,vsCol,drop=F]
  mDest
}
