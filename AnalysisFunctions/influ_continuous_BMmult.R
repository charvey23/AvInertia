## Custom version of sensiPhy::influ_continuous()
## For use with multivariate data to calculate multivariate
## BM rates via geomorph::compare.evol.rates()
##
## Modified on 2021-04-23 from influ_continuous() from sensiPhy package
## version 0.8.5 to use compare.evol.rates() from geomorph 4.0.0 by V.B.B.
##
## The only model option for this function is "BM"
## "sigsq" values in $sensi.estimates will all be multivariate sigma sq values


influ_continuous_BMmult <-
  function(data,
           phy,
           model = "BM",
           gp,
           iter,
           method = "simulation",
           cutoff = 2,
           track = TRUE,
           ...) {
    {
      ## Model must be BM
      model <- "BM"

      if (is.null(model))
        stop("model must be specified, e.g. 'OU' or 'lambda'")
      full.data <- data
      phy <- phy
      N <- dim(full.data)[1]
      mod.0 <-
        geomorph::compare.evol.rates(
          A = full.data,
          phy = phy,
          method = method,
          gp = gp,
          iter = 999,
          ...
        )

      sigsq.0 <- mod.0$sigma.d.all
      if (model == "BM") {
        optpar.0 <- NA
      }
      sensi.estimates <-
        data.frame(
          species = numeric(),
          sigsq = numeric(),
          DIFsigsq = numeric(),
          sigsq.perc = numeric()
        )
      counter <- 1
      errors <- NULL
      if (track == TRUE)
        pb <- utils::txtProgressBar(min = 0,
                                    max = N,
                                    style = 3)
      for (i in 1:N) {
        crop.data <- full.data[c(1:N)[-i], ]
        crop.phy <- ape::drop.tip(phy, setdiff(phy$tip.label,
                                               rownames(crop.data)))
        crop.gp <- gp[-i]
        mod = try(geomorph::compare.evol.rates(
          A = crop.data,
          phy = crop.phy,
          method = method,
          gp = crop.gp,
          iter = 999,
          ...
        )
        ,
        TRUE)

        if (isTRUE(class(mod) == "try-error")) {
          error <- i
          names(error) <- rownames(full.data$data)[i]
          errors <- c(errors, error)
          next
        }
        else {
          sp <- phy$tip.label[i]
          sigsq <- mod$sigma.d.all
          DIFsigsq <- sigsq - sigsq.0
          sigsq.perc <- round((abs(DIFsigsq / sigsq.0)) * 100,
                              digits = 1)
          if (model == "BM") {
            optpar <- NA
          }
          DIFoptpar <- optpar - optpar.0
          optpar.perc <- round((abs(DIFoptpar / optpar.0)) *
                                 100, digits = 1)
          if (track == TRUE)
            utils::setTxtProgressBar(pb, i)
          estim.simu <- data.frame(
            sp,
            sigsq,
            DIFsigsq,
            sigsq.perc,
            stringsAsFactors = F
          )
          sensi.estimates[counter,] <- estim.simu
          counter = counter + 1
        }
      }
      if (track == TRUE)
        on.exit(close(pb))
      sDIFsigsq <-
        sensi.estimates$DIFsigsq / stats::sd(sensi.estimates$DIFsigsq)
      sensi.estimates$sDIFsigsq <- sDIFsigsq
      if (model == "BM") {
        sDIFoptpar <- NA
      }
      param0 <- list(sigsq = sigsq.0, optpar = optpar.0)
      reorder.on.sigsq <-
        sensi.estimates[order(abs(sensi.estimates$sDIFsigsq),
                              decreasing = T), c("species", "sDIFsigsq")]
      influ.sp.sigsq <-
        as.character(reorder.on.sigsq$species[abs(reorder.on.sigsq$sDIFsigsq) >
                                                cutoff])
      if (model == "BM") {
        influ.sp.optpar <-
          "No optpar calculated for BM-model. Influential species not calculated"
      }
      res <-
        list(
          call = match.call(),
          cutoff = cutoff,
          data = full.data,
          optpar = model,
          full.model.estimates = param0,
          influential.species = list(influ.sp.sigsq = influ.sp.sigsq,
                                     influ.sp.optpar = influ.sp.optpar),
          sensi.estimates = sensi.estimates,
          errors = errors
        )
      class(res) <- "sensiInflu.TraitEvol"
      if (length(res$errors) > 0) {
        warning("Some species deletion presented errors, please check: output$errors")
      }
      else {
        res$errors <- "No errors found."
      }
      return(res)
    }
  }
