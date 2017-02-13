

bin_plot <- function(data, binwidth, l_units){
  newDat <- bin_mat(data, binwidth)
  
  LB_pars <- new("LB_pars", default = FALSE)
  # LB_pars@Species <- "SPECIES"
  LB_obj <- new("LB_obj")
  
  LB_lengths <- new("LB_lengths", 
                    file = newDat,
                    dataType = "freq")
  
  LB_lengths@Years <- as.numeric(gsub("X", "", colnames(newDat)[-c(1:2)]))
  LB_lengths@LData <- as.matrix(newDat[, -c(1, 2)])
  LB_lengths@LMids <- newDat$lmidp
  LB_lengths@L_units <- l_units
  plotSize(LB_lengths)
}

bin_mat <- function(data, binwidth) {
  # First column is the current length class
  # Remaining columns are years
  # Returns data frame:
  # lclass, lmidp, and years
  
  current_binwidth <- data[2,1] - data[1,1] 
  
  if(current_binwidth > binwidth) {
    stop("Bin width (", binwidth, ") should be greater than original bin width (",
         current_binwidth, ").")
  }
  
  data <- as.data.frame(data,
                        StringsAsFactors = FALSE)
  
  data[is.na(data)] <- 0
  
  minCL <- floor((min(data[,1]) - .5) / binwidth) * binwidth
  maxCL <- ceiling((max(data[,1]) + .5) / binwidth) * binwidth
  
  minCL <- ifelse(minCL <= 0, 0, minCL)
  
  break_list <- seq(minCL,
                    maxCL, 
                    binwidth)
  
  data$LC <- cut(data[,1], 
                 breaks = break_list,
                 include.lowest = T)
  
  dWide <- aggregate(data[, 3:ncol(data)-1], 
                     by = list(data$LC),
                     sum, na.rm = FALSE)
  
  dWide <- merge(dWide, 
                 data.frame(LC = as.character(levels(data$LC))),
                 by.x = "Group.1",
                 by.y = "LC",
                 all = TRUE,
                 sort = FALSE)

  dWide[is.na(dWide)] <- 0
  
  ints <- seq(minCL + binwidth / 2,
              maxCL + binwidth / 2,
              binwidth)
  
  dWide <- data.frame(lclass = dWide[, 1],
                      lmidp = ints[-length(ints)],
                      dWide[, 2:ncol(dWide)])
  return(dWide) 
}

lb_ind <- function(data, 
                   binwidth,
                   linf, 
                   lmat,
                   weight) {
  
  # if(is.null(weight)) {
  #   message(paste0("Note: without weight, Lmaxy and associated reference points cannot be calculated.",
  #                  " The plots will still show up."))
  # } else {
  weight <- bin_mat(weight, binwidth)
  # }  
  newDat <- bin_mat(data, binwidth)
  
  cols <- as.numeric(gsub("X", "", colnames(newDat)[-c(1:2)]))
  startyear <- min(cols)
  endyear <- max(cols)
    
  ind_names <- c("Year", "L75", "L25", "Lmed",
                 "L90", "L95", "Lmean", "Lc",
                 "LFeM", "Lmaxy", "Lmat", "Lopt",
                 "Linf", "Lmax5", "Lmean_LFeM",
                 "Lc_Lmat", "L25_Lmat", "Lmean_Lmat",
                 "Lmean_Lopt", "L95_Linf", 
                 "Lmaxy_Lopt", "Lmax5_Linf", "Pmega", "Pmegaref")
  Ind <- data.frame(matrix(ncol = length(ind_names), 
                           nrow = endyear - startyear + 1))
  names(Ind) <- ind_names
  Ind$Year <- startyear:endyear
  
  #  regrouping with selected length class width
  longDat <- melt(newDat[, -1],
                  id.var = "lmidp",
                  value.name = "number",
                  variable.name = "year")
  longDat$year <- as.numeric(gsub("X", "", as.character(longDat$year)))

  Year <- seq(startyear, endyear)
  res <- data.frame(year = Year,
                    lmidp = NA,
                    nmax = NA, 
                    lc = NA)
  
  # newDat <- bin_mat(data, binwidth)

  for (j in 3:ncol(newDat)) {
    for (i in 2:nrow(newDat)) {
      if(newDat[i + 1, j] - newDat[i, j] >= 0) {
        next
      } else {
        res$lmidp[j-2] <- newDat$lmidp[i]
        res$nmax[j-2] <- newDat[i,j]
        a <- res$nmax[j - 2] / 2
        df1 <- newDat[, c(2, j)]
        for (k in 1:nrow(df1)) {
          if (df1[k, 2] < a) {
            next
          } else {  
            res$lc[j - 2] <- df1[k, 1]
          }
          break
        }
      }
      break
    }
  }
  
  Ind$Lc <- res$lc
  Ind$Lmat <- lmat
  Ind$Lopt <- 2/3 * linf
  Ind$Linf <- linf
  
  final <- newDat[,-1]
  for(jj in (1:length(Year)) + 1){
    j <- jj-1 
    
    final2 <- final[, c(1, jj)]
    colnames(final2) <- c("lngth","number")
    
    final2$cumsum <- cumsum(final2[, 2])
    final2$cumsum_perc <- final2$cumsum / sum(final2$number)
    
    # find mean top 5% 
    # from largest starting
    numb <- as.data.frame(final2[rev(order(final2$lngth)), "number"])
    colnames(numb) <- "number"
    numb$cum <- cumsum(numb$number) 
    numb$lngth <- final2[rev(order(final2$lngth)),"lngth"] 
    numb$cumperc <- round(numb$cum/sum(numb$number),5)  
    numb$num5 <- 0
    numb[numb$cumperc <= 0.05, "num5"] <- numb[numb$cumperc <= 0.05, "number"]
    numb[max(which(numb$cumperc <= 0.05)) + 1, 
         "num5"] <- (0.05 - numb[max(which(numb$cumperc <= 0.05)),
                                 "cumperc"]) * sum(numb$number)
    Ind[j, "Lmax5"] <- sum(numb$num5 * numb$lngth) / sum(numb$num5)
    
    # indicators
    Ind[j, "L75"] <- min(final2[which(final2$cumsum_perc >= 0.75), "lngth"])
    Ind[j, "L25"] <- min(final2[which(final2$cumsum_perc >= 0.25), "lngth"])
    Ind[j, "Lmed"] <- min(final2[which(final2$cumsum_perc >= 0.5), "lngth"])
    Ind[j, "L95"] <- min(final2[which(final2$cumsum_perc >= 0.95), "lngth"])
    Ind[j, "L90"] <- min(final2[which(final2$cumsum_perc >= 0.90), "lngth"])

    # calculate mean of individuals above Lc
    final3 <- final2[final2$lngth >= Ind[j, "Lc"], ]
    Ind[j, "Lmean"] <- sum(final3$lngth * final3$number) / sum(final3$number)
    
    # length class with max yield
    if(!is.null(weight)) {
    final2$biomass <- final2$number * weight[, jj]
    Ind[j, "Lmaxy"] <- final2[final2$biomass == max(final2$biomass), "lngth"]
    } else {
    Ind[j, "Lmaxy"] <- NA
    }
    
    Lopt <- (2 / 3) * linf
    
    # proportion larger Lopt+10%
    Ind[j, "Pmega"] <- sum(final2[which(final2$lngth >= (Lopt + 0.1 * Lopt)),
                                  "number"]) / sum(final2$number)
    Ind[j, "Year"] <- Year[j]
    Ind[j, "Pmegaref"] <- 0.3   # proxy reference point of 30% in catch
    Ind[j, "LFeM"] <- 0.75 * Ind[j, "Lc"] + 0.25 * Ind[j, "Linf"]
  }
  
  #calculate various ratios
  Ind$Lmaxy_Lopt <- Ind$Lmaxy / Ind$Lopt
  Ind$L95_Linf <- Ind$L95 / Ind$Linf
  Ind$Lmean_LFeM <- Ind$Lmean / Ind$LFeM
  Ind$Lmean_Lmat <- Ind$Lmean / Ind$Lmat
  Ind$Lmean_Lopt <- Ind$Lmean / Ind$Lopt
  Ind$Lmax5_Linf <- Ind$Lmax5 / Ind$Linf
  Ind$Lc_Lmat <- Ind$Lc / Ind$Lmat
  Ind$L25_Lmat <- Ind$L25 / Ind$Lmat
  return(Ind)
}

lb_plot <- function(data,
                    binwidth,
                    l_units,
                    linf,
                    lmat,
                    weight) {

  Ind <- lb_ind(data = data,
                binwidth = binwidth,
                linf = linf,
                lmat = lmat,
                weight = weight)

  ymax_a <- max(Ind$L95,
                Ind$Lmax5,
                Ind$Lmean,
                Ind$Lc,
                Ind$Linf,
                Ind$L25,
                na.rm = TRUE)
  ymax_b <- max(Ind$Lmax5_Linf,
                Ind$L95_Linf,
                Ind$Pmega,
                Ind$Pmegaref,
                Ind$Lc_Lmat,
                Ind$L25_Lmat,
                na.rm = TRUE)
  ymax_c <- max(Ind$L75,
                Ind$Lmean,
                Ind$Lopt,
                Ind$Lmaxy,
                Ind$Lmat,
                Ind$L25,
                na.rm = TRUE)
  ymax_d <- max(Ind$Lmean_Lopt,
                Ind$Lmaxy_Lopt,
                na.rm = TRUE)
  ymax_e <- max(Ind$Lmat,
                Ind$Lmean,
                Ind$LFeM,
                na.rm = TRUE)
  ymax_f <- max(Ind$Lmean_LFeM,
                na.rm = TRUE)

    par(mar = c(4, 4, 3, 6.1),
        mfrow = c(3, 2),
        family = "serif",
        cex = 1,
        xpd = TRUE)
    
  par(mfg = c(1,1))
  plot(Linf ~ Year,
       data = Ind,
       ylab = paste("Length (", l_units, ")", sep = ""),
       col = "transparent",
       main = "(a) Conservation",
       xlab = "Year",
       ylim = c(0,
                ymax_a * 1.1),
       bty = "l")
  legend("topright",
         legend = c(expression(L["95%"]),
                    expression(L["max5%"]),
                    expression(L["mat"]),
                    expression(L["c"]),
                    expression(L["inf"]),
                    expression(L["25%"])),
         lwd = c(2,2,1,2,1,1),
         lty = c(1,1,3,1,3,1),
         text.col = c("purple", #L95 solid
                      "black", #Lmax5 Solid
                      "grey40", #Lmat dash
                      "blue", #Lc solid
                      "black", #Linf dash
                      "red"), # L25 solid
         col = c("purple",
                 "black",
                 "grey40",
                 "blue",
                 "black",
                 "red"),
         bty = "n",
         seg.len = 0.7,
         inset = c(-0.15, 0))
  lines(L95 ~ Year,
        data = Ind,
        lwd = 2,
        col = "purple")
  lines(Lmax5 ~ Year,
        data = Ind,
        lwd = 2,
        col = "black")
  lines(Lmat ~ Year,
        data = Ind,
        lwd = 1,
        col = "grey40",
        lty = "dashed")
  lines(Lc ~ Year,
        data = Ind,
        lwd = 2,
        col = "blue")
  lines(Linf ~ Year,
        data = Ind,
        lwd = 1,
        col = "black",
        lty = "dashed")
  lines(L25 ~ Year,
        data = Ind,
        lwd = 1,
        col = "red")

  par(mfg = c(2, 1))
  plot(Linf ~ Year,
       data = Ind,
       ylab = paste("Length (", l_units, ")", sep = ""),
       main = "(c) Optimal Yield",
       col = "transparent",
       xlab = "Year",
       ylim = c(0,
                ymax_c * 1.1),
       bty="l")
  lines(L75 ~ Year,
        data = Ind,
        lwd = 1,
        col = "red")
  lines(Lmean ~ Year,
        data = Ind,
        lwd = 2,
        col = "darkred")
  lines(Lopt ~ Year,
        data = Ind,
        lwd = 1,
        col = "black",
        lty = "dashed")
  lines(Lmaxy ~ Year,
        data = Ind,
        lwd = 2,
        col = "green")
  lines(Lmat ~ Year,
        data = Ind,
        lwd = 1,
        col = "grey40",
        lty = "dashed")
  lines(L25 ~ Year,
        data = Ind,
        lwd = 1,
        col = "red")
  legend(x = "topright",
         legend = c(expression(L["75%"]),
                    expression(L["mean"]),
                    expression(L["opt"]),
                    expression(L["maxy"]),
                    expression(L["mat"]),
                    expression(L["25%"])),
         text.col = c("red",
                      "darkred",
                      "black",
                      "green",
                      "grey40",
                      "red"),
         col = c("red",
                 "darkred",
                 "black",
                 "green",
                 "grey40",
                 "red"),
         lty = c(1, 1, 3, 1, 3, 1),
         lwd = c(1, 2, 1, 2, 1, 1),
         bty = "n",
         seg.len = 0.7,
         inset = c(-0.15, 0))

  par(mfg = c(3,1))
  plot(Lmat ~ Year,
       data = Ind,
       type = "n",
       ylab = paste("Length (", l_units, ")", sep = ""),
       main = "(e) Maximum Sustainable Yield",
       xlab = "Year",
       ylim = c(0,
                ymax_e * 1.1),
       bty = "l")
  lines(Lmat ~ Year,
        data = Ind,
        lwd = 1,
        col = "grey40",
        lty = "dashed")
  lines(Lmean ~ Year,
        data = Ind,
        lwd = 2,
        col = "darkred")
  lines(LFeM ~ Year,
        data = Ind,
        lwd = 2,
        col = "blue",
        lty = "dashed")
  legend(x = "topright",
         legend = c(expression(L["mat"]),
                    expression(L["mean"]),
                    expression(L["F=M"])),
         text.col = c("grey40",
                      "darkred",
                      "blue"),
         col = c("grey40",
                 "darkred",
                 "blue"),
         lty = c(3, 1, 3),
         lwd = c(1, 2, 2),
         bty = "n",
         seg.len = 0.7,
         inset = c(-0.15, 0))

  par(mfg = c(1, 2))
  plot(Lmax5_Linf ~ Year,
       data = Ind,
       type = "n",
       ylab = "Indicator Ratio",
       col = "transparent",
       main = "(b) Conservation",
       xlab = "Year",
       ylim = c(0,
                ymax_b * 1.1),
       bty = "l")
  lines(Lmax5_Linf ~ Year,
        data = Ind,
        lwd = 2,
        col = "black")
  lines(L95_Linf ~ Year,
        data = Ind,
        lwd = 2,
        col = "purple")
  lines(Pmega ~ Year,
        data = Ind,
        lwd = 2,
        col = "blue")
  lines(Pmegaref ~ Year,
        data = Ind,
        lwd = 1,
        col = "black",
        lty="dashed")
  lines(Lc_Lmat ~ Year,
        data = Ind,
        lwd = 2,
        col = "red")
  lines(L25_Lmat ~ Year,
        data = Ind,
        lwd = 2,
        col = "darkred")
  legend("topright",
         legend = c(expression(L["max5%"]/L["inf"]),
                    expression(L["95%"]/L["inf"]),
                    expression(P["mega"]),
                    expression("30%"),
                    expression(L["c"]/L["mat"]),
                    expression(L["25"]/L["mat"])),
         lwd = c(2, 2, 2, 1, 2, 2),
         lty = c(1, 1, 1, 3, 1, 1),
         text.col = c("black",
                      "purple",
                      "blue",
                      "black",
                      "red",
                      "darkred"),
         col = c("black",
                 "purple",
                 "blue",
                 "black",
                 "red",
                 "darkred"),
         bty = "n",
         seg.len = 0.7,
         inset = c(-0.2, 0))

  par(mfg = c(2,2))
  plot(Lmean_Lopt ~ Year,
       data = Ind,
       ylab = "Indicator Ratio",
       col = "transparent",
       main = "(d) Optimal yield",
       xlab = "Year",
       type = "n",
       ylim = c(0, ymax_d * 1.1),
       bty = "l")
  lines(Lmean_Lopt ~ Year,
        data = Ind,
        lwd = 2,
        col = "darkred")
  lines(Lmaxy_Lopt ~ Year,
        data = Ind,
        lwd = 2,
        col = "green")
  legend("topright",
         legend = c(expression(L["mean"]/L["opt"]),
                    expression(L["maxy"]/L["opt"])),
         col = c("darkred",
                 "green"),
         text.col = c("darkred",
                      "green"),
         lty = c(1, 1),
         lwd = c(1, 1),
         bty = "n",
         seg.len = 0.7,
         inset = c(-0.2, 0))

  par(mfg = c(3,2))
  plot(Lmean_LFeM ~ Year,
       data = Ind,
       type = "n",
       ylab = "Indicator Ratio",
       col = "transparent",
       main = "(f) Maximum sustainable yield",
       xlab = "Year",
       ylim = c(0,
                ymax_f * 1.1),
       bty = "l")
  lines(Lmean_LFeM ~ Year,
        data = Ind,
        lwd = 2,
        col = "blue")
  legend("topright",
         legend = expression(L["mean"]/L["F=M"]),
         col = "blue",
         text.col = "blue",
         lty = 1,
         lwd = 1,
         bty = "n",
         seg.len = 0.7,
         inset = c(-0.2, 0))
# dev.off()
}

lb_doc <- function(data,
                   binwidth,
                   l_units,
                   linf,
                   lmat,
                   stock,
                   weight,
                   filename) {
  base_text_prop <- textProperties(font.family = "Calibri")


  tmpfile_lfd <-  tempfile(pattern = paste(stock ,"_LFD_", sep = ""),
                           fileext = ".png")
  tmpfile_indicator <-  tempfile(pattern = paste(stock ,"_Indicator_", sep = ""), 
                                 fileext = ".png")
  

  ggsave(plot = suppressMessages(bin_plot(data, 
                                          binwidth,
                                          l_units)), 
         filename = tmpfile_lfd,
         width = 16,
         height = 15, 
         units = "cm",
         dpi = 400)
  
  png(tmpfile_indicator,
      bg = "white",
      pointsize = 5,
      units = "cm",
      width = 16,
      height = 15,
      res = 400)
  
  lb_plot(data,
          binwidth,
          l_units,
          linf,
          lmat,
          weight)
  dev.off()
  
  subscript_pot <- function(text) {
    pot(text, 
        textProperties(font.family = "Calibri", 
                       font.weight = "normal", 
                       font.size = 11,
                       vertical.align = "subscript"))
  }

  head_pot <- function(text) {
    pot(text, 
        textProperties(font.family = "Calibri", 
                       font.weight = "normal", 
                       font.size = 11))
  }
  
  linf_pot <- head_pot("L") + subscript_pot("inf") + head_pot(paste0(" = ", linf))
  lmat_pot <- head_pot("L") + subscript_pot("mat") + head_pot(paste0(" = ", lmat))
  
  draftDoc <- docx(template = "data/report_template.docx",
                   title = stock)

  draftDoc <- addParagraph(draftDoc, 
                           value = paste(stock, "Length Based Indicators", sep = " "),
                           text.properties = base_text_prop,
                           bookmark = "STOCK_NAME")
  
  draftDoc <- addParagraph(draftDoc, 
                           value = lmat_pot,
                           text.properties = base_text_prop,
                           bookmark = "LMAT")
  
  draftDoc <- addParagraph(draftDoc, 
                           value = linf_pot,
                           text.properties = base_text_prop,
                           bookmark = "LINF")
  
  
  draftDoc <- addFlexTable(draftDoc,
                           flextable = lb_table(data,
                                                binwidth,
                                                l_units,
                                                linf,
                                                lmat,
                                                weight),
                           bookmark = "INDICATOR_TABLE")
  
  draftDoc <- deleteBookmark(draftDoc,
                             bookmark = "INDICATOR_TABLE")
  
  draftDoc <- addImage(draftDoc,
                       filename = tmpfile_indicator,
                       width = 17/2.54,
                       height = 16/2.54,
                       bookmark = "INDICATOR_PLOT")
  
  draftDoc <- addImage(draftDoc, 
                       filename = tmpfile_lfd,
                       width = 12/2.54,
                       height = 10/2.54, 
                       bookmark = "LFD_PLOT")
  
  draftDoc <- addFlexTable(draftDoc,
                           flextable = lb_raw_dat(data = data,
                                                binwidth,
                                                l_units),
                           bookmark = "LFD_TABLE")

  draftDoc <- deleteBookmark(draftDoc,
                             bookmark = "LFD_TABLE")
  

    draftDoc <- addFlexTable(draftDoc,
                           flextable = lb_raw_dat(data = weight,
                                                binwidth,
                                                l_units),
                           bookmark = "WAL_TABLE")
  
  draftDoc <- deleteBookmark(draftDoc,
                             bookmark = "WAL_TABLE")

  
  writeDoc(draftDoc, file = filename)
  file.remove(tmpfile_indicator)
  file.remove(tmpfile_lfd)
}

lb_table <- function(data,
                     binwidth,
                     l_units,
                     linf,
                     lmat,
                     weight) {
  
  
  Ind <- lb_ind(data = data,
                binwidth = binwidth,
                linf = linf,
                lmat = lmat,
                weight = weight)
  
  
  ref_level <- c(0, 1, 1, 0.8, 0.3, 0.9, 1)
  years <- seq(max(Ind$Year) - 2,  max(Ind$Year))
  flex_dat <- Ind[c("Year",
                    "Lc_Lmat", "L25_Lmat", "Lmax5_Linf", "Pmega",
                    "Lmean_Lopt", "Lmean_LFeM")]
  flex_dat <- flex_dat[flex_dat$Year %in% years,]
  flex_dat <- round(flex_dat, 2)
  
  flex_tab <- FlexTable(flex_dat, header.columns = FALSE,
                        body.text.props = textProperties(font.family = "Calibri", 
                                                         font.weight = "normal", 
                                                         font.size = 9))
  
  flex_tab <- addHeaderRow(flex_tab, 
                           text.properties = textNormal(),
                           value = c("",
                                     "Conservation",
                                     "Optimizing Yield",
                                     "MSY"), 
                           colspan = c(1, 4, 1, 1))
  
  flex_tab <- addHeaderRow(flex_tab, 
                           text.properties = textNormal(),
                           value = c("Year",
                                     rep("", 6)), 
                           colspan = rep(1,7))
  
  head_pot <- function(text) {
    pot(text, 
        textProperties(font.family = "Calibri", 
                       font.weight = "normal", 
                       font.size = 9))
  }
  
  subscript_pot <- function(text) {
    pot(text, 
        textProperties(font.family = "Calibri", 
                       font.weight = "normal", 
                       font.size = 9,
                       vertical.align = "subscript"))
  }
  
  flex_tab[1:2, , to = "header"] <- textProperties(font.family = "Calibri", 
                                                   font.weight = "normal", 
                                                   font.size = 9)
  flex_tab[2, 2, to = "header"] <- head_pot("L") + subscript_pot("c") + head_pot(" / L") + subscript_pot("mat")
  flex_tab[2, 3, to = "header"] <- head_pot("L") + subscript_pot("25%") + head_pot(" / L") + subscript_pot("mat")
  flex_tab[2, 4, to = "header"] <- head_pot("L") + subscript_pot("max 5") + head_pot(" / L") + subscript_pot("inf")
  flex_tab[2, 5, to = "header"] <- head_pot("P") + subscript_pot("mega")
  flex_tab[2, 6, to = "header"] <- head_pot("L") + subscript_pot("mean") + head_pot(" / L") + subscript_pot("opt")
  flex_tab[2, 7, to = "header"] <- head_pot("L") + subscript_pot("mean") + head_pot(" / L") + subscript_pot("F = M")
  flex_tab[1:2, , to = "header"] <- parProperties(text.align = "center",
                                                  padding = 3)
  flex_tab[1:2, , to = "header"] <- cellProperties(background.color = "#E8EAEA")
  
  # Body formatting
  flex_tab[ , 1 , to = "body"] <- parProperties(text.align = "center",
                                                padding = 3 )
  flex_tab[ , 2:7 , to = "body"] <- parProperties(text.align = "right",
                                                  padding = 3 )
  
  flex_tab <- setFlexTableBackgroundColors(flex_tab, 
                                           j = 2,
                                           colors = ifelse(flex_dat[,2] > ref_level[2],
                                                           "#aec640",
                                                           "#f15d2a"))
  flex_tab <- setFlexTableBackgroundColors(flex_tab, 
                                           j = 3,
                                           colors = ifelse(flex_dat[,3] > ref_level[3],
                                                           "#aec640",
                                                           "#f15d2a"))
  flex_tab <- setFlexTableBackgroundColors(flex_tab, 
                                           j = 4,
                                           colors = ifelse(flex_dat[,4] > ref_level[4],
                                                           "#aec640",
                                                           "#f15d2a"))
  flex_tab <- setFlexTableBackgroundColors(flex_tab, 
                                           j = 5,
                                           colors = ifelse(flex_dat[,5] > ref_level[5],
                                                           "#aec640",
                                                           "#f15d2a"))
  flex_tab <- setFlexTableBackgroundColors(flex_tab, 
                                           j = 6,
                                           colors = ifelse(flex_dat[,6] > ref_level[6],
                                                           "#aec640",
                                                           "#f15d2a"))
  flex_tab <- setFlexTableBackgroundColors(flex_tab, 
                                           j = 7,
                                           colors = ifelse(flex_dat[,7] >= ref_level[7],
                                                           "#aec640",
                                                           "#f15d2a"))
  setFlexTableWidths(flex_tab, c((1.4/2.54), rep(2.76/2.54, 6)))
  
  return(flex_tab)
}

lb_raw_dat <- function(data,
                       binwidth,
                       l_units,
                       weight) {
  
  Ind <- bin_mat(data = data,
                 binwidth = binwidth)
  
  Ind[,-1] <- round(Ind[,-1], 2)
  cols <- c(paste0("Length class \n(", l_units, ")"),
             paste0("Length midpoint \n(", l_units, ")"),
             gsub("X", "", colnames(Ind[-(1:2)])))
  
  flex_tab <- FlexTable(Ind, header.columns = FALSE,
                        body.text.props = textProperties(font.family = "Calibri", 
                                                         font.weight = "normal", 
                                                         font.size = 9),
                        body.par.props = parProperties(text.align = "right",
                                                       padding = 3))
  
  flex_tab <- addHeaderRow(flex_tab, 
                           text.properties = textProperties(font.family = "Calibri", 
                                                            font.weight = "normal", 
                                                            font.size = 9),
                           par.properties = parProperties(text.align = "center",
                                                          padding = 3),
                           value = cols)
  
  # Body formatting
  flex_tab[ , 1:2 , to = "body"] <- parProperties(text.align = "center",
                                                  padding = 3 )
  return(flex_tab)
}
