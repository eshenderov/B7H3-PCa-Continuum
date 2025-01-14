
df_transpose <- function(df){
  col_names <- rownames(df)
  row_names <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  df <- as.data.frame(df, row.names = row_names)
  colnames(df) <- col_names
  return(df)
}

transform_raw_counts <- function(obj, assay, slot, transformation, output_slot) {


  return(obj)

}
