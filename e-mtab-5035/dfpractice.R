# creating data frame
data_frame <- data.frame(col1 = rep(letters[1:4], each = 2),
                         col2 = 1:8
)
print("Original DataFrame")
print(data_frame)

# assigning row names to data frame
rownames(data_frame) <- c("row1","row2","row3","row4",
                          "row5", "row6","row7","row8")

# getting rows
rows <- c("row1","row3","row5","row8")

# extracting data frame rows
data_mod <- data_frame[rownames(data_frame) %in% rows, ]
print("Modified DataFrame")
print(data_mod)
