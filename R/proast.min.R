# Concatenated files for minimal proast package
#' Run dose-response modeling using PROAST.
#'
#' @param interactive_mode A TRUE/FALSE value specifying whether you want to run interactively (i.e., TRUE, the default) or using command-line mode (i.e., FALSE, non-interactive)
#' @param datatype Options: continuous, individual data.
#' @return Results from PROAST.
f.proast <- function(odt = list(), ans.all = 0, er = FALSE, resize = FALSE, 
    scale.ans = FALSE, const.var = F, show.warnings = F, interactive_mode = TRUE,
    datatype = NULL, model_choice = NULL, setting_choice = NULL, nested_model_choice = NULL,
    indep_var_choice = NULL, Vyans_input = NULL, covariates = NULL, custom_CES = 0.05, model_selection = NULL, lower_dd = NULL,
    upper_dd = NULL, selected_model = NULL, adjust_CES_to_group_SD = NULL, model_averaging = NULL, num_bootstraps = NULL, display_plots = TRUE) {
    message(paste0("Independent variable: ", indep_var_choice))
    if (interactive_mode == TRUE) {
        message("Running in interactive mode...")
    } else {
        message("Running in batch/command line mode...")
        message(paste0("Independent variable: ", indep_var_choice))
        results_env <- new.env()
    }
    if (exists("track")) 
        print("f.proast")
    if (is.list(odt)) 
        if (length(odt$EXP) > 0 || length(odt$MLE) > 0) {
            cat("you probably forgot the comma in f.proast()\n")
            ans.cont <- menu(c("no", "yes"), title = "Do you want to continue nonetheless?")
            if (ans.cont == 1) 
                return(invisible())
        }
    odt.name <- deparse(substitute(odt))
    f.overlap()
    quit <- F
    if (mode(ans.all) == "list") {
        ans.all <- f.adjust.saved(ans.all)
        ans.all$const.var <- const.var
        ans.all$quick.ans <- 1
        if (ans.all$cont) 
            quit <- f.con(ans.all, list.logic = T, indep_var_choice = indep_var_choice, Vyans_input = Vyans_input, covariates = covariates, custom_CES = custom_CES, model_selection = model_selection, lower_dd = lower_dd, upper_dd = upper_dd, interactive_mode = interactive_mode, selected_model = selected_model, adjust_CES_to_group_SD = adjust_CES_to_group_SD, model_averaging = model_averaging, num_bootstraps = num_bootstraps, results_env = results_env, display_plots = display_plots)
        else quit <- f.cat(ans.all, list.logic = T)
        if (quit) {
            if (interactive_mode == FALSE) {
                result <- as.list(results_env)
                result_df <- parse_PROAST_output(result)
                print(result_df)
                return(list(result, result_df))
            } else {
                cat("\n\n   Goodbye,\n    PROAST \n\n")
                return(invisible())
            }
        }
    }
    message(paste0("Independent variable: ", indep_var_choice))
    if (!er) {
        ans.all <- f.ini(odt)
        ans.all$odt.name <- odt.name
        ans.all$scale.ans <- scale.ans
        ans.all$const.var <- const.var
        ans.all$show.warnings <- show.warnings
        f.graphwin.size(resize)
        if (interactive_mode == TRUE) {
            dtype <- menu(c("continuous, individual data", "binary", 
            "ordinal", "quantal", "clustered continuous, individual data", 
            "clustered quantal", "continuous, summary data", 
            "clustered continuous, summary data", "quantal, CxT", 
            "other"), title = "\nWhat type of response data do you want to consider?")
        } else {
            dtype_options <- setNames(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                         c("continuous, individual data", "binary", "ordinal", 
                           "quantal", "clustered continuous, individual data", 
                           "clustered quantal", "continuous, summary data", 
                           "clustered continuous, summary data", "quantal, CxT", 
                           "other"))
            dtype <- dtype_options[datatype]
            message(paste0("Setting data type to ", names(dtype_options[datatype])))
        }
        if (dtype == 10) 
            dtype <- eval(parse(prompt = "Give value for dtype > "))
        if (dtype == 7) 
            dtype <- 10
        if (dtype == 8) 
            dtype <- 15
        if (dtype == 9) 
            dtype <- 84
        ans.all$dtype <- dtype
        if (dtype %in% c(1, 5, 15, 25, 26)) 
            ans.all$plt.mns <- 1
        if (dtype %in% c(10, 15, 250, 260)) {
            ans.all$plt.mns <- 3
            ans.all$CI.plt <- T
        }
        f.assign(".Pr.last", ans.all)

        # Add mappings for both dtype and quick.ans menus
        model_options <- c("single model", "change settings first", 
        "select model 3 or 5 from various families of models", 
        "select model 3 from various nested families of models", 
        "select model 5 from various nested families of models", 
        "select model 15 in terms of RPF")

        nested_model_options <- c("single model", "", 
        "select model 3 or 5 from various nested families of models", 
        "select model 3 from various families of models", 
        "select model 5 from various families of models", 
        "select model 15 in terms of RPF")
        
        setting_options <- c("single model", "set of models")

        if (dtype %in% c(1, 5, 10, 15)) {
            ans.all$cont <- TRUE
            if (interactive_mode == TRUE) {
                ans.all$quick.ans <- menu(model_options, title = "\nDo you want to fit a single model or fit various nested families of models?")
            } else {
                if (!is.null(model_choice)) {
                    ans.all$quick.ans <- which(model_options == model_choice)
                    message(paste0("Selected model option: '", model_choice, "' (", ans.all$quick.ans, ")"))
                } else {
                    stop("Parameter 'model_choice' must be provided in non-interactive mode.")
                }
            }
            if (ans.all$quick.ans == 2) {
                ans.all$change <- rep(F, ans.all$nrQ)
                message(paste0("Independent variable: ", indep_var_choice))
                message(ans.all)
                ans.all <- f.change.settings(ans.all, choose = T, indep_var_choice = indep_var_choice)
                message(ans.all)
                if (interactive_mode == TRUE) {
                    ans.all$quick.ans <- menu(nested_model_options, title = "\nDo you want to fit a single model or fit various nested families of models?")
                } else {
                    if (!is.null(nested_model_choice)) {
                        ans.all$quick.ans <- which(nested_model_options == nested_model_choice)
                        message(paste0("Selected nested model option: '", nested_model_choice, "' (", ans.all$quick.ans, ")"))
                    } else {
                        stop("Parameter 'nested_model_choice' must be provided in non-interactive mode.")
                    }
                }
            }
            if (ans.all$quick.ans > 1) 
                quit <- f.con(ans.all, indep_var_choice = indep_var_choice, Vyans_input = Vyans_input, covariates = covariates, custom_CES = custom_CES, model_selection = model_selection, lower_dd = lower_dd, upper_dd = upper_dd, interactive_mode = interactive_mode, selected_model = selected_model, adjust_CES_to_group_SD = adjust_CES_to_group_SD, model_averaging = model_averaging, num_bootstraps = num_bootstraps, results_env = results_env, display_plots = display_plots)
        }
        if (dtype %in% c(2, 3, 4, 6, 84)) {
            ans.all$cont <- FALSE
            if (interactive_mode == TRUE) {
                ans.all$quick.ans <- menu(c("single model", "change settings first", "set of models"),
                                          title = "\nDo you want to fit a set of models, or choose a single model?")
                if (ans.all$quick.ans == 2) {
                    ans.all$change <- rep(F, ans.all$nrQ)
                    message(paste0("Independent variable: ", indep_var_choice))
                    message(ans.all)
                    ans.all <- f.change.settings(ans.all, choose = T, indep_var_choice = indep_var_choice)
                    ans.all$quick.ans <- menu(c("single model", "set of models"), 
                    title = "\nDo you want to fit a set of models, or choose a single model?")
                    if (ans.all$quick.ans == 2) {
                    ans.all$quick.ans <- 3
                    ans.all <- f.quick.cat(ans.all, interactive_mode = interactive_mode, results_env = results_env)
                    }
                }
            } else {
                if (!is.null(setting_choice)) {
                    ans.all$quick.ans <- which(setting_options == setting_choice)
                    message(paste0("Selected model option: '", setting_choice, "' (", ans.all$quick.ans, ")"))
                } else {
                    stop("Parameter 'setting_choice' must be provided in non-interactive mode.")
                }
            }
            if (ans.all$quick.ans == 3) 
                quit <- f.cat(ans.all)
            }
        if (quit) {
            if (interactive_mode == FALSE) {
                result <- as.list(results_env)
                result_df <- parse_PROAST_output(result)
                print(result_df)
                return(list(result, result_df))
            } else {
                cat("\n\n   Goodbye,\n    PROAST \n\n")
                return(invisible())
            }
        }
        ans.all$change[1] <- T
        ans.all <- f.change.settings(ans.all, indep_var_choice = indep_var_choice)
        f.assign(".Pr.last", ans.all)
        x <- ans.all$data.0[, ans.all$xans]
        x <- x[!is.na(x)]
        if (is.list(x)) 
            for (jj in 1:length(x)) f.check.nonneg.num(x[, jj])
        else f.check.nonneg.num(x)
        ans.all$change[2] <- T
        ans.all <- f.change.settings(ans.all, indep_var_choice = indep_var_choice)
        if (!is.null(ans.all$data.0)) 
            ans.all$odt.tmp <- ans.all$data.0
        else ans.all$odt.tmp <- ans.all$odt
        ans.all$data.0 <- f.firstcheck(ans.all)
        f.assign(".Pr.last", ans.all)
    }
    if (er) {
        ans.all <- .Pr.last
        if (ans.all$xans[1] == 0 || ans.all$yans[1] == 0) {
            cat("\nOption er=T cannot be used, start analysis all over\n\n")
            return(invisible())
        }
        ans.all$quick.ans <- 1
    }
    if (ans.all$quick.ans == 1) {
        ans.all <- f.execute(ans.all)
    }
    if (ans.all$cont) 
        quit <- f.con(ans.all, indep_var_choice = indep_var_choice, Vyans_input = Vyans_input, covariates = covariates, CES = CES, model_selection = model_selection, lower_dd = lower_dd, upper_dd = upper_dd, interactive_mode = interactive_mode, selected_model = selected_model, adjust_CES_to_group_SD = adjust_CES_to_group_SD, model_averaging = model_averaging, num_bootstraps = num_bootstraps, results_env = results_env, display_plots = display_plots)
    else quit <- f.cat(ans.all)
    if (quit) {
        if (interactive_mode == FALSE) {
            result <- as.list(results_env)
            result_df <- parse_PROAST_output(result)
            print(result_df)
            return(list(result, result_df))
        } else {
            cat("\n\n   Goodbye,\n    PROAST \n\n")
            return(invisible())
        }
    }
}



f.overlap <- function() {
    if (exists("track")) 
        print("f.overlap")
    attached <- search()
    attached <- substring(attached, 9, 14)
    lst1 <- attached == "proast"
    attached <- substring(search(), 1, 4)
    lst2 <- attached == "file"
    pos.pr <- 0
    if (any(lst1)) {
        pos.pr <- (1:length(lst1))[lst1]
    }
    if (any(lst2)) {
        pos.pr <- (1:length(lst2))[lst2]
    }
    pos.pr <- pos.pr[1]
    if (pos.pr > 0) {
        if (exists("last.fit", where = pos.pr)) {
            rm(last.fit, pos = pos.pr)
            if (exists("newfit")) 
                rm(newfit, pos = pos.pr)
        }
        if (length(pos.pr) > 1) 
            cat("\n\nATTENTION: you have attached more than one PROAST package\n\n")
        obj.pr <- objects(pos = pos.pr)
        obj.pr[obj.pr == "ans.all"] <- ""
        wdir <- getwd()
        dum <- substring(wdir, 30, 36)
        if (dum != "proast") {
            obj.wd <- objects(pos = 1)
            obj.wd <- obj.wd[obj.wd != "newfit"]
            lst <- obj.wd %in% obj.pr
            overlap <- obj.wd[lst]
            if (sum(lst) > 0) {
                cat("\nProast cannot be run !")
                cat("\npress escape to stop this session\n")
                cat("\nand then remove or rename the following object(s):\n")
                for (ii in 1:sum(lst)) {
                  cat(overlap[ii])
                  cat("\n")
                }
                f.press.key.to.continue(NA)
                stop()
            }
        }
    }
    if (exists("track")) 
        print("f.overlap:   END")
    return(invisible())
}



f.ini <- function(odt = NULL, gui = FALSE, no.seed = FALSE) {
    if (exists("track")) 
        print("f.ini")
    ans.all <- list()
    if (!no.seed) 
        ans.all$seed.bt <- 125
    else ans.all$seed.bt <- 0
    ans.all$PRversion <- f.version()
    nrQ <- 27
    ans.all$gui <- gui
    if (is.list(odt) && length(odt$nvar) == 1) {
        ans.all$varnames <- odt$varnames
        ans.all$nvar <- odt$nvar
        ans.all$odt <- odt$data
        ans.all$info <- odt$info
    }
    else {
        ans.all$varnames <- colnames(odt)
        ans.all$nvar <- length(ans.all$varnames)
        ans.all$odt <- odt
        ans.all$info <- " "
    }
    ans.all$data.0 <- ans.all$odt
    ans.all$nrQ <- nrQ
    ans.all$change <- rep(F, nrQ + 1)
    ans.all$xans <- 0
    ans.all$yans <- 0
    ans.all$sf.x <- 1
    ans.all$sf.y <- 1
    ans.all$nans <- 0
    ans.all$CED.model <- F
    ans.all$sans <- 0
    ans.all$sd.se <- 0
    ans.all$detlim <- 0
    ans.all$Vdetlim <- 0
    ans.all$ans.detlim <- 1
    ans.all$detlim.col <- 0
    ans.all$auto.detlim <- F
    ans.all$transf.ans <- 3
    ans.all$cens.up <- NA
    ans.all$nest.no <- 0
    ans.all$fct4.no <- 0
    ans.all$fct5.no <- 0
    ans.all$covar.dd <- 0
    ans.all$regr.resid <- NA
    ans.all$pred.value <- NA
    ans.all$x.mn <- NA
    ans.all$ans.m6.sd <- 1
    ans.all$sign.q <- 1
    ans.all$const.var <- FALSE
    ans.all$opposing <- 1
    ans.all$nr.models <- 2
    ans.all$out.ans <- 2
    ans.all$first.outliers <- TRUE
    ans.all$model.switch <- 0
    ans.all$tans <- 0
    ans.all$sf.t <- 1
    ans.all$model.type <- 1
    ans.all$P.gof <- 0.05
    ans.all$aic.crit <- 2
    ans.all$aic <- NA
    ans.all$ces.ans <- 3
    ans.all$alfa.length <- 0
    ans.all$alfa.start <- 10
    ans.all$alfa.Pvalue <- NA
    ans.all$nr.gr <- 1
    ans.all$CI.plt <- F
    ans.all$shift.tmp <- 0
    ans.all$cc.inf <- 1e+18
    ans.all$logit.ma <- FALSE
    ans.all$remove.no <- 0
    ans.all$remove.lev <- list()
    ans.all$remove.name <- "none"
    ans.all$select.no <- 0
    ans.all$select.lev <- list()
    ans.all$select.name <- "all"
    ans.all$fct1 <- 1
    ans.all$fct2 <- 1
    ans.all$fct3 <- 1
    ans.all$fct4 <- 1
    ans.all$fct5 <- 1
    ans.all$fct1.no <- 0
    ans.all$fct2.no <- 0
    ans.all$fct3.no <- 0
    ans.all$fct4.no <- 0
    ans.all$fct5.no <- 0
    ans.all$displ.no <- 0
    ans.all$covar.no <- 0
    ans.all$nr.aa <- 0
    ans.all$nr.bb <- 0
    ans.all$nr.var <- 0
    ans.all$nr.cc <- 0
    ans.all$nr.dd <- 0
    ans.all$incr.decr.no <- 0
    ans.all$ref.lev <- NA
    ans.all$model.ans <- 0
    ans.all$MLE <- NA
    ans.all$converged <- 0
    ans.all$lb <- NA
    ans.all$ub <- NA
    ans.all$scale.ans <- FALSE
    ans.all$CES <- 0
    ans.all$NES.ans <- 1
    ans.all$BMR <- 0
    ans.all$CED <- NA
    ans.all$m6.ED50 <- FALSE
    ans.all$x1 <- 0
    ans.all$x2 <- 0
    ans.all$CED.ref <- NA
    ans.all$BMD <- NA
    ans.all$CED.matr <- matrix(NA, ncol = 2, nrow = 1)
    ans.all$ES.abs <- 2
    ans.all$plot.type <- 2
    ans.all$l.ty <- 1
    ans.all$heading <- ""
    ans.all$x.leg <- ""
    ans.all$y.leg <- ""
    ans.all$t.leg <- ""
    ans.all$low.y <- NA
    ans.all$upp.y <- NA
    ans.all$x.lim <- rep(NA, 2)
    ans.all$y.lim <- NA
    ans.all$shift <- 0
    ans.all$min.x <- NA
    ans.all$max.x <- NA
    ans.all$plt.mns <- 1
    ans.all$cex.1 <- NA
    ans.all$cex.2 <- NA
    ans.all$identify <- FALSE
    ans.all$label <- 0
    ans.all$ans.plt <- 0
    ans.all$down <- F
    ans.all$CI.dum <- 1
    ans.all$LVM <- T
    ans.all$logprob <- F
    ans.all$do.MA <- F
    ans.all$MA.running <- F
    ans.all$combi.ans <- F
    ans.all$svg.plots <- F
    ans.all$ans.scale <- 0
    ans.all$bounds.ans <- 2
    ans.all$modelname <- ""
    assign(".ypos", 0, immediate = T, pos = 1)
    ans.all$store.name <- 0
    ans.all$list.logic <- F
    ans.all$CES.cat <- 1
    ans.all$def.exc <- NA
    ans.all$ced.lines <- F
    ans.all$categ.ans <- 0
    ans.all$cat.ans <- 0
    ans.all$reverse.ans <- 0
    ans.all$trace <- T
    ans.all$group <- 0
    ans.all$nruns <- 0
    ans.all$constr.steepness <- 0.01
    ans.all$no.CI <- F
    ans.all$interrupt <- TRUE
    if (gui) 
        ans.all$interrupt <- FALSE
    ans.all$rank.low <- 0
    ans.all$cc.fix <- NA
    ans.all$increase <- -1
    ans.all$increase.00 <- -1
    ans.all$vv <- 1
    ans.all$hill <- 0
    ans.all$WAPP <- FALSE
    ans.all$efsa.tool <- FALSE
    ans.all$simu <- FALSE
    ans.all$output <- TRUE
    ans.all$name.wapp <- NA
    ans.all$plotprefix <- NA
    color <- 1:100
    color <- color[!color %in% c(7, 15)]
    ans.all$color <- color
    ans.all$mark <- c(2, 4:25, 32:255)
    ans.all$x <- NA
    ans.all$y <- NA
    ans.all$yy <- NA
    ans.all$yans.denom <- 0
    ans.all$pooled.sd <- NA
    ans.all$mn.log <- 0
    ans.all$sd2.log <- 0
    ans.all$y.add <- 0
    ans.all$nr.dd <- 1
    ans.all$nr.cc <- 1
    ans.all$kk <- NA
    ans.all$ttt <- 0
    ans.all$ttt.txt <- ""
    ans.all$lump.ans <- 2
    ans.all$scores.orig <- 0
    ans.all$displ.fact <- 1
    ans.all$factor.name <- character()
    ans.all$covariate <- 1
    ans.all$fct1.txt <- ""
    ans.all$fct2.txt <- ""
    ans.all$fct3.txt <- ""
    ans.all$fct4.txt <- ""
    ans.all$fct5.txt <- ""
    ans.all$gr.txt <- ""
    ans.all$covar.txt <- ""
    ans.all$fct1.name <- ""
    ans.all$fct2.name <- ""
    ans.all$fct3.name <- ""
    ans.all$fct4.name <- ""
    ans.all$fct5.name <- ""
    ans.all$covar.name <- ""
    ans.all$twice <- F
    ans.all$par.start <- NA
    ans.all$eps <- 1e-12
    ans.all$fit.ans <- 1
    ans.all$fitted <- F
    ans.all$new.model <- F
    ans.all$scale.dum <- 0
    ans.all$text.par <- character()
    ans.all$text.par.0 <- NA
    ans.all$nrp <- 0
    ans.all$npar <- NA
    ans.all$adjust.start <- F
    ans.all$cond.ans <- 3
    ans.all$lst.control <- f.control(3)
    ans.all$loglik <- NA
    ans.all$low.yln <- NA
    ans.all$show <- ""
    ans.all$nl <- "\n"
    ans.all$tb <- "\t"
    ans.all$sp <- " "
    ans.all$sep <- FALSE
    ans.all$par.boot <- 0
    ans.all$boot <- F
    ans.all$xy.lim <- NA
    ans.all$y.lim.CI <- NA
    ans.all$pi.full <- NA
    ans.all$x.full <- NA
    ans.all$fct1.full <- NA
    ans.all$fct2.full <- NA
    ans.all$kk.tot <- 0
    ans.all$nn.tot <- 0
    ans.all$regr.start <- numeric()
    ans.all$full6.done <- FALSE
    ans.all$th.start <- NA
    ans.all$th.0.start <- 0
    ans.all$sig.start <- 0.25
    ans.all$regr.par <- numeric()
    ans.all$regr.par.matr <- matrix()
    ans.all$regr.par.matr.bb <- matrix()
    ans.all$th.par <- 0
    ans.all$th.par.vec <- NA
    ans.all$nth <- 1
    ans.all$sig.par <- NA
    ans.all$decr.zz <- T
    ans.all$response.matr <- 0
    ans.all$plotcurve <- FALSE
    conf.int <- matrix(NA, ncol = 2)
    ans.all$conf.int <- conf.int
    ans.all$conf.lev <- 0.9
    ans.all$quick.ans <- 1
    if (gui) 
        ans.all$quick.ans <- 3
    ans.all$ab.txt <- ""
    ans.all$npar.full <- NA
    ans.all$npar.null <- NA
    ans.all$loglik.full <- NA
    ans.all$loglik.null <- NA
    ans.all$aic.full <- NA
    ans.all$aic.null <- NA
    ans.all$full.ans <- 1
    if (ans.all$gui) {
        ans.all$data <- data.frame()
        ans.all$data.is.from.file <- FALSE
        ans.all$last.upload.dir <- ""
        ans.all$data.name <- ""
        ans.all$comma.for.decimal <- f.use.comma.for.decimal()
    }
    ans.all$summ.lst <- list()
    ans.all$notes <- "ATTENTION"
    ans.all$report.pars <- NA
    ans.all$alert.full <- NA
    ans.all$show.no.trend <- NA
    ans.all$DA.ans <- 1
    class(ans.all) <- c("list", "proast.session")
    if (exists("track")) 
        print("f.ini : END")
    return(ans.all)
}



f.graphwin.size <- function(resize = F) {
    if (exists("track")) 
        print("f.graphwin.size")
    if (!exists(".gw.size")) {
        size.x <- 5
        size.y <- 5
        f.assign(".gw.size", c(size.x, size.y))
    }
    if (resize) {
        f.create.graphwin(.gw.size[1], .gw.size[2])
        oke <- F
        while (!oke) {
            gw.x <- .gw.size[[1]][1]
            gw.y <- .gw.size[[2]][1]
            cat("\ncurrent values for horizontal and vertical size of the graphical window are:\n")
            cat(gw.x, " and ", gw.y, "\n")
            gw.x <- eval(parse(prompt = paste("\nGive value for horizontal size > ")))
            gw.y <- eval(parse(prompt = paste("Give value for vertical size > ")))
            f.create.graphwin(gw.x, gw.y)
            oke <- 1 - menu(c("no", "yes"), title = "\nIs this window size okay?")
            lst <- dev.list()
            curr <- dev.cur()
            lst <- lst[lst != curr]
            for (ii in lst) dev.off(ii)
        }
        assign(".gw.size", c(gw.x, gw.y), pos = 1, immediate = T)
    }
}



f.assign <- function(name, obj) {
    if (exists("track")) 
        print("f.assign")
    obj$date <- date()
    if (!is.character(name)) 
        print("f.assign:  names and object are reversed")
    assign(name, obj, pos = 1, immediate = T)
    if (exists("track")) 
        print("f.assign:   END")
}



f.change.settings <- function(ans.all, choose = F, indep_var_choice = NULL, interactive_mode = TRUE) {
  message(indep_var_choice)
    if (exists("track")) 
        print("f.change.settings")
    ans.all$choose <- choose
    with(ans.all, {
        if (0) {
            if (quick.ans == 1 && change[2]) {
                if (dtype %in% c(10, 15, 250, 260)) 
                  if (sans[1] == 0) 
                    change[3] <- T
                if (dtype %in% c(4, 6, 84)) 
                  if (nans[1] == 0) 
                    change[8] <- T
            }
        }
        while (choose || sum(change[1:(nrQ - 1)]) != 0) {
            if (choose) {
                settings <- f.define.settings(quick.ans, dtype, 
                  xans, nrQ, model.ans)
                if ((cont || model.type == 2) && model.ans %in% 
                  c(32, 45)) 
                  settings[c(14, 16, 17)] <- ""
                set.ans <- menu(settings, title = "\nWhich setting(s) do you want to change?")
                if (!(set.ans %in% c(9, 18, 23, nrQ))) 
                  ans.all$fitted <- F
                change[set.ans] <- T
            }
            if (change[1]) {
              message("Independent variable:")
              message(indep_var_choice)
              message("moving on...")
              if (!is.null(indep_var_choice)) {
                message(paste0("Using ", indep_var_choice, " as the independent variable."))
                if (indep_var_choice %in% ans.all$varnames) {
                    ans.all$xans <- which(ans.all$varnames == indep_var_choice)
                } else {
                    stop(paste("The specified independent variable choice '", indep_var_choice, "' is not one of the column names."))
                }
                change[1] <- F  # Set change[1] to FALSE after successful setting of independent variable.
              } else {
                  cat("\nQ1: Which variable do you want to consider as the independent variable?\n    (e.g. dose, age)\n")
                  cat(" ")
                  cat(paste(1:nvar, ":", varnames, "\n"))
                  xans <- eval(parse(prompt = paste(" -------- > ")))
                  ans.all$xans <- xans
                  change[1] <- F
                  DA.ans <- 2
                  if (length(xans) > 1) {
                    ans.all$displ.no <- menu(varnames[1:nvar], 
                      title = "\nQ9: Which column defines the single doses and mixtures?\n    (type 0 if none)")
                    title.tmp <- "Do you want to fit dose addition model, or parallel curves model? \n         ATTENTION: the parallel curve model cannot be chosen when the data include mixtures"
                    DA.ans <- menu(c("dose addition model", "parallel curves model"), 
                      title = title.tmp)
                    ans.all$gr.txt <- levels(factor(data.0[, ans.all$displ.no]))
                    if (DA.ans == 2) {
                      ans.all <- f.add.dosecols(ans.all)
                      ans.all$DA.ans <- 2
                      ans.all$fct2.no <- ans.all$displ.no
                    }
                }
              }
            }
            if (change[2]) {
                yans <- menu(varnames[1:nvar], title = "\nQ2: Which response variable do you want to analyse ? \n")
                ans.all$yans <- yans
                change[2] <- F
                if (dtype %in% c(4, 6, 84)) 
                  change[8] <- T
                if (dtype %in% c(10, 15, 250, 260)) 
                  change[3] <- T
                if (dtype == 5) {
                  ans.all$nest.no <- menu(varnames[1:nvar], title = "\nGive column number of the nested factor \n                         (type 0 to skip nested analysis)")
                }
                ans.all$detlim <- 0
                ans.all$Vdetlim <- 0
                ans.all$detlim.col <- 0
                ans.all$reverse.ans <- 0
            }
            if (change[3]) {
                if (dtype == 3) {
                  cat.ans <- menu(1:nth, title = paste("\nQ3:  give severity category for dichotomization \n"))
                  ans.all$cat.ans <- cat.ans
                }
                else {
                  sans <- menu(varnames[1:nvar], title = paste("\nQ3a:  number of the variation statistic associated to", 
                    varnames[yans], " \n"))
                  sd.se <- menu(c("standard deviations", "standard errors"), 
                    title = "\nQ3b: Do you have standard deviations or standard errors associated with the means?")
                  nans <- menu(varnames[1:nvar], title = "\nQ3c: Give the associated sample size\n")
                  ans.all$sans <- sans
                  ans.all$sd.se <- sd.se
                  ans.all$nans <- nans
                }
                change[3] <- F
            }
            if (change[4]) {
                go.on <- 1
                count <- 0
                while (go.on > 0) {
                  remove.no <- menu(varnames[1:nvar], title = paste("\nQ4: If you want to remove outliers, type number of factor to indicate them", 
                    "\n--- type 0 if you do not want to remove any (more) data --- \n"))
                  if (remove.no == 0) {
                    go.on <- 0
                    tmp.lev <- "none"
                  }
                  else {
                    count <- count + 1
                    remove.fact <- data.0[, remove.no]
                    if (!is.numeric(remove.fact)) {
                      levels.all <- sort(levels(factor(remove.fact)))
                      nr.lev <- length(levels.all)
                      cat("\nQ4b: Give number(s) before the level(s) you want to remove\n")
                      cat(" ")
                      cat(paste(1:nr.lev, ":", levels.all, "\n"))
                      tmp.no <- eval(parse(prompt = paste(" -------- > ")))
                      tmp.lev <- levels.all[tmp.no]
                    }
                    if (is.numeric(remove.fact)) {
                      levels.all <- unique(remove.fact)
                      cat("\n the values of this variable are: \n")
                      print(sort(levels.all))
                      cat("\nQ4b: Give values of the levels you want to remove\n")
                      tmp.lev <- eval((parse(prompt = paste("   > "))))
                      tmp.lev <- levels.all[levels.all %in% tmp.lev]
                    }
                    ans.all$remove.lev <- sort(tmp.lev)
                    ans.all$remove.no <- remove.no
                    ans.all$data.0 <- data.0
                    data.0 <- f.subset(ans.all, remove = T)
                    ans.all$data.0 <- data.0
                    remove.tmp <- paste(paste(varnames[remove.no], 
                      ":", paste(tmp.lev, collapse = " ")), sep = "\n")
                    if (remove.name == "none") 
                      remove.name <- character()
                    remove.name <- paste(remove.name, "\n  ", 
                      remove.tmp)
                  }
                  if (count == 0) {
                    cat("\nATTENTION: complete dataset will be restored\n")
                    data.0 <- odt
                    remove.name <- character()
                    select.name <- character()
                    if (length(unique(Vdetlim) == 1)) {
                      Vdetlim <- rep(Vdetlim[1], length(data.0[, 
                        1]))
                      data.0 <- cbind(data.0, Vdetlim)
                    }
                    ans.all$data.0 <- data.0
                  }
                  if (length(remove.name) == 0) 
                    remove.name <- "none"
                  ans.all$remove.name <- remove.name
                }
                ans.all$covar.no <- 0
                change[4] <- F
            }
            if (change[5]) {
                yans.denom <- menu(varnames[1:nvar], title = "\nWhich variable do you want as denominator ? (0 = NONE)\n")
                change[5] <- F
                ans.all$yans.denom <- yans.denom
            }
            if (change[6]) {
                if (yans != 0) {
                  y <- data.0[, yans]
                  print(sort(y))
                  cat("\nthe lowest nonzero value is: ", min(y[y > 
                    0]), "\n")
                  if (min(y) > 0) 
                    cat("\n\nNOTE: lowest observation is positive \n\n")
                }
                detlim <- eval(parse(prompt = paste("\nQ6: Give global value of detection limit, \nor type 0 for individual values > ")))
                if (detlim == 0) {
                  ans.all$detlim.col <- menu(varnames[1:nvar], 
                    title = "\nGive column with individual detection limits\n")
                  ans.all$Vdetlim <- data.0[, ans.all$detlim.col]
                }
                if (detlim.col == 0) {
                  ans.all$Vdetlim <- rep(detlim, length(x))
                  ans.all$data.0 <- cbind(data.0, Vdetlim)
                  ans.all$detlim.col <- nvar + 1
                  ans.all$nvar <- nvar + 1
                  ans.all$varnames[ans.all$nvar] <- "detlim"
                }
                change[6] <- F
            }
            if (change[7]) {
                y <- data.0[, yans]
                y <- y[!is.na(y)]
                cat("\nthe lowest nonzero observation is: ", 
                  min(y[y > 0]), "\n")
                y.add <- eval(parse(prompt = "\n\nQ7: Give small value to be added to the data > "))
                change[7] <- F
                ans.all$y.add <- y.add
            }
            if (change[8]) {
                ans.all$nans <- menu(varnames[1:nvar], title = "\nQ3b: Give the associated sample sizes")
                if (dtype == 84) 
                  ans.all$tans <- menu(varnames[1:nvar], title = "\nQ3c: Give number of time variable")
                change[8] <- F
            }
            if (change[9]) {
                ans.all$displ.no <- menu(varnames[1:nvar], title = "\nQ9: Which variable do you want to consider for distinct plotting?\n    (type 0 if none)")
                change[9] <- F
            }
            if (change[10]) {
                xy.ans <- menu(c("dose", "response"), title = "\nQ10: Do you want to scale dose or response? \n")
                if (xy.ans == 1) {
                  if (ans.all$xans != 0) {
                    x <- data.0[, xans]
                    cat("\n(Nonzero) doses range from ", min(x[x != 
                      0], na.rm = T), "to ", max(x))
                  }
                  ans.all$sf.x <- eval(parse(prompt = "\nQ10a: Give scaling factor for dose > "))
                }
                if (xy.ans == 2) 
                  ans.all$sf.y <- eval(parse(prompt = "\nQ10b: Give scaling factor for response > "))
                ans.all <- f.clear(ans.all)
                if (tans > 0) {
                  ttt <- data.0[, tans]
                  cat("\n\n\n(Nonzero) times range from ", min(ttt[ttt != 
                    0]), "to ", max(ttt))
                  sf.t <- eval(parse(prompt = "\nQ10: Give scaling factor for time > "))
                  ans.all$sf.t <- sf.t
                }
                change[10] <- F
                if (model.ans == 46) {
                  cat("\nscaling factor on dose not implemented for model 46, it will be ignored\n")
                  f.press.key.to.continue()
                }
            }
            if (change[11]) {
                go.on <- 1
                count <- 0
                while (go.on > 0) {
                  select.no <- menu(varnames[1:nvar], title = "\nQ11a: Give number of (another) factor for which you want to select data\n                          --- type 0 if you want to analyze all (remaining) data --- \n")
                  if (select.no == 0) {
                    go.on <- 0
                    tmp.lev <- "none"
                  }
                  else {
                    count <- count + 1
                    select.fact <- data.0[, select.no]
                    if (!is.numeric(select.fact)) {
                      levels.all <- sort(levels(factor(select.fact)))
                    }
                    if (is.numeric(select.fact)) {
                      levels.all <- sort(unique(select.fact))
                    }
                    nr.lev <- length(levels.all)
                    cat("\nQ11b: Give number(s) associated with the level(s) you want to select\n")
                    cat(" ")
                    cat(paste(1:nr.lev, ":", levels.all, "\n"))
                    tmp.no <- eval(parse(prompt = paste(" -------- > ")))
                    tmp.lev <- levels.all[tmp.no]
                    zz.1 <- data.0[, select.no]
                    zz.2 <- zz.1[!is.na(zz.1)]
                    if (quick.ans == 1 && length(zz.1) > length(zz.2)) 
                      cat("\n\nATTENTION: there are missing values in this factor\n")
                    ans.all$select.no <- select.no
                    ans.all$select.lev <- tmp.lev
                    ans.all$data.0 <- data.0
                    data.0 <- f.subset(ans.all, select = T)
                    ans.all$data.0 <- data.0
                    name.tmp <- paste(varnames[select.no], paste(tmp.lev, 
                      collapse = " "))
                    if (select.name == "all") 
                      select.name <- character()
                    select.name <- paste(select.name, "\n  ", 
                      name.tmp)
                  }
                }
                if (count == 0) {
                  cat("\nATTENTION: complete dataset will be restored\n")
                  data.0 <- odt
                  select.name <- character()
                  remove.name <- character()
                  if (length(unique(Vdetlim) == 1)) {
                    Vdetlim <- rep(Vdetlim[1], length(data.0[, 
                      1]))
                    data.0 <- cbind(data.0, Vdetlim)
                  }
                  ans.all$data.0 <- data.0
                }
                if (length(select.name) == 0) 
                  select.name <- "all"
                change[11] <- F
                ans.all$select.lev <- tmp.lev
                ans.all$select.name <- select.name
                ans.all$covar.no <- 0
            }
            if (change[12]) {
                if (!cont) {
                  if (dtype == 2) 
                    lump.ans <- menu(c("yes", "no"), title = "\nQ12: Do you want to consider observed fractions per dose group\n \n                  i.e. you have at least 5 observations per dose group?")
                  else lump.ans <- 2
                  ans.all$lump.ans <- lump.ans
                  ans.all$yans <- 0
                }
                change[12] <- F
            }
            if (change[13]) {
                fct1.no <- menu(varnames[1:nvar], title = "\nQ13: Give number of factor serving as covariate with respect to parameter a \n --- type 0 if none --- \n  ")
                change[13] <- F
                ans.all$fct1.no <- fct1.no
                ans.all$fct1.name <- varnames[fct1.no]
            }
            if (change[14]) {
                fct2.no <- menu(varnames[1:nvar], title = "\nQ14: Give number of factor serving as covariate with respect to parameter b \n --- type 0 if none --- \n  ")
                change[14] <- F
                ans.all$fct2.no <- fct2.no
                ans.all$fct2.name <- varnames[fct2.no]
            }
            if (change[15]) {
                if (dtype == 4) 
                  fct3.no <- menu(varnames[1:nvar], title = "\nQ15: Give number of factor serving as covariate with respect to parameter theta \n --- type 0 if none --- \n")
                else if (dtype == 3) 
                  fct3.no <- menu(varnames[1:nvar], title = "\nQ15: Give number of factor serving as covariate with respect to    \n                    parameter var \n --- type 0 if none --- \n")
                else if (dtype == 6) 
                  fct3.no <- menu(varnames[1:nvar], title = "\nQ15: Give number of factor serving as covariate with respect to    \n                    parameter alpha \n --- type 0 if none --- \n")
                else fct3.no <- menu(varnames[1:nvar], title = "\nQ15: Give number of factor serving as covariate with respect to    \n                    parameter var \n --- type 0 if none --- \n")
                change[15] <- F
                ans.all$fct3.no <- fct3.no
                ans.all$fct3.name <- varnames[fct3.no]
            }
            if (change[16]) {
                if (model.ans %in% c(2, 4, 7, 9, 12, 14, 17, 
                  19, 22, 24)) {
                  cat("\nATTENTION: parameter d does not occur in this model\n\n")
                  f.press.key.to.continue()
                }
                if (model.type == 2 && model.ans %in% c(13, 15, 
                  23, 25)) {
                  cat("\nATTENTION: covariate on parameter d not implemented for LVMs in terms of CED \n\n")
                  f.press.key.to.continue()
                  return(ans.all)
                }
                if (model.ans != 57) 
                  fct5.no <- menu(varnames[1:nvar], title = "\nQ16: Give number of factor serving as covariate with respect to\n    parameter d \n --- type 0 if none --- \n  ")
                else fct5.no <- menu(varnames[1:nvar], title = "\nQ16: Give number of factor serving as covariate with respect to\n    parameter b1 \n --- type 0 if none --- \n  ")
                change[16] <- F
                ans.all$fct5.no <- fct5.no
                ans.all$fct5.name <- varnames[fct5.no]
            }
            if (change[17]) {
                if (model.ans %in% c(2, 3, 7, 8, 12, 13, 17, 
                  18, 22, 23)) {
                  cat("\nATTENTION: parameter c does not occur in this model\n\n")
                  f.press.key.to.continue()
                }
                if (model.type == 2 && model.ans %in% c(15, 25)) {
                  cat("\nATTENTION: covariate on parameter c not implemented for LVMs in terms of CED \n\n")
                  f.press.key.to.continue()
                  return(ans.all)
                }
                fct4.no <- menu(varnames[1:nvar], title = "\nQ16: Give number of factor serving as covariate with respect to\n    parameter c \n --- type 0 if none --- \n  ")
                change[17] <- F
                ans.all$fct4.no <- fct4.no
                ans.all$fct4.name <- varnames[fct4.no]
            }
            if (change[18]) {
                if (is.na(x[1])) 
                  ans.all <- f.execute(ans.all)
                if (dtype %in% c(1, 5, 10, 15, 25, 26)) {
                  res.MS <- with(ans.all, res.MS <- f.con.tst(ans.all), 
                    return(res.MS))
                  ans.all$res.MS <- res.MS
                }
                else if (dtype == 4) {
                  f.cat.tst(ans.all)
                }
                else cat("\n\nNOAEL approach not implemented for dtype", 
                  dtype, "\n\n")
                change[18] <- F
            }
            if (change[19] && cont) 
                if (dtype %in% c(1, 5, 10, 15, 25, 26, 250, 260)) {
                  transf.ans <- menu(c("no transformation", "sqrt-transformation", 
                    "log-transformation"), title = "\nQ19: Which transformation do you want?")
                  ans.all$transf.ans <- transf.ans
                  if (dtype %in% c(10, 15, 250, 260)) 
                    switch(transf.ans, dtype <- 250, dtype <- 260, 
                      dtype <- 10)
                  if (dtype %in% c(1, 25, 26)) 
                    switch(transf.ans, dtype <- 25, dtype <- 26, 
                      dtype <- 1)
                  if (dtype == 5) 
                    switch(transf.ans, dtype <- 25, dtype <- 26, 
                      dtype <- 5)
                  if (dtype %in% c(1, 5, 10, 15)) 
                    cat("\nData will be analysed on log-scale\n")
                  if (dtype %in% c(25, 250)) 
                    cat("\nData will be analysed on original scale\n")
                  if (dtype %in% c(26, 260)) 
                    cat("\nData will be analysed on square root-scale\n")
                  ans.all$dtype <- dtype
                  if (ans.all$yans > 0) 
                    ans.all <- f.execute(ans.all)
                  change[19] <- F
                }
                else cat("\nATTENTION : data cannot be transformed\n")
            if (change[19] && !cont) {
                ans.all$constr.steepness <- eval(parse(prompt = paste("Q19: Give value for lower constraint on the steepness parameter > ")))
                change[19] <- F
            }
            if (change[20]) {
                ans.all$aic.crit <- eval(parse(prompt = paste("Q20: Give value for \n      critical difference in AICs among models (default = 2) > ")))
                change[20] <- F
            }
            if (change[21]) {
                if (cont) {
                  if (yans != 0) {
                    y <- data.0[, yans]
                    print(sort(y))
                  }
                  cens.up <- eval(parse(prompt = paste("\n\nQ21: Give value for right censoring > ")))
                }
                if (dtype %in% c(4, 6)) {
                  cat(paste(1:nvar, ":", varnames, "\n"))
                  cens.no <- eval(parse(prompt = paste("\n\nQ21: Give column number indicating the censoring > ")))
                  cens.up <- data.0[, cens.no]
                }
                change[21] <- F
                ans.all$cens.up <- cens.up
                change[21] <- F
            }
            if (change[22]) {
                ans.all$fit.ans <- menu(c("error in response only (default)", 
                  "error in both response and dose (minimize sum of products)"), 
                  title = "\nQ22: What type of regression do you want?")
                change[22] <- F
            }
            if (change[23]) {
                ans.all$cond.ans <- menu(c("mild conditions", 
                  "moderate conditions", "strict conditions (R default)"), 
                  title = "\nQ23: How strict do you want the fitting conditions to be?")
                ans.all$lst.control <- f.control(ans.all$cond.ans)
                change[23] <- F
                change[nrQ] <- T
            }
            if (change[24]) {
                ans.all$tans <- menu(varnames[1:nvar], title = "\nQ24: Give column number of time variable \n")
                change[24] <- F
            }
            if (change[25]) {
                ans.all$conf.lev <- eval(parse(prompt = paste(" give confidence level --- > ")))
                change[25] <- F
            }
            if (change[26]) {
                if (model.ans %in% c(0, 6, 15, 25, 52, 54, 46)) 
                  ans.all$incr.decr.no <- menu(varnames[1:nvar], 
                    title = "\nQ24: Give column number indicating increasing or decreasing dose-response \n")
                else cat("\n\n this option is not possibe for this model, choose model 6, 15, 25, 52, 54 or 46 \n\n")
                change[26] <- F
            }
            if (change[27]) {
                choose <- F
            }
            ans.all$show <- ""
            ans.all$change <- change
            ans.all$full6.done <- FALSE
            f.assign(".Pr.last", ans.all)
        }
        if (exists("track")) 
            print("f.change.settings:  END")
        return(ans.all)
    })
}



f.remove.NAs <- function(xans = 0, tans = 0, yans = 0, sans = 0, nans = 0, covar.no = 0, 
    fct1.no = 0, fct2.no = 0, fct3.no = 0, fct4.no = 0, fct5.no = 0, 
    dfr, output = TRUE) {
    if (exists("track")) 
        print("f.remove.NAs")
    if (length(xans) == 1) 
        lst.na.x <- is.na(dfr[, xans])
    else {
        Mx <- dfr[, xans]
        nr.dosecol <- length(xans)
        Mx <- !is.na(Mx)
        Mx.sum <- rowSums(Mx)
        lst.na.x <- Mx.sum != nr.dosecol
    }
    lst.na.t <- rep(F, length(dfr[, 1]))
    lst.na.y <- rep(F, length(dfr[, 1]))
    lst.na.sd <- rep(F, length(dfr[, 1]))
    lst.na.nn <- rep(F, length(dfr[, 1]))
    lst.na.covar <- rep(F, length(dfr[, 1]))
    lst.na.fct1 <- rep(F, length(dfr[, 1]))
    lst.na.fct2 <- rep(F, length(dfr[, 1]))
    lst.na.fct3 <- rep(F, length(dfr[, 1]))
    lst.na.fct4 <- rep(F, length(dfr[, 1]))
    lst.na.fct5 <- rep(F, length(dfr[, 1]))
    if (tans != 0) 
        lst.na.t <- is.na(dfr[, tans])
    if (yans != 0) 
        lst.na.y <- is.na(dfr[, yans])
    if (sans != 0) 
        lst.na.sd <- is.na(dfr[, sans])
    if (nans != 0) 
        lst.na.nn <- is.na(dfr[, nans])
    if (covar.no != 0) 
        lst.na.covar <- is.na(dfr[, covar.no])
    if (fct1.no != 0) 
        lst.na.fct1 <- is.na(dfr[, fct1.no])
    if (fct2.no != 0) 
        lst.na.fct2 <- is.na(dfr[, fct2.no])
    if (fct3.no != 0) 
        lst.na.fct3 <- is.na(dfr[, fct3.no])
    if (fct4.no != 0) 
        lst.na.fct4 <- is.na(dfr[, fct4.no])
    if (fct5.no != 0) 
        lst.na.fct5 <- is.na(dfr[, fct5.no])
    if (output) {
        if (sum(lst.na.x) > 0) 
            cat("\n\nThere are ", sum(lst.na.x), "missing values in the x column\n")
        if (sum(lst.na.t) > 0) 
            cat("\n\nThere are ", sum(lst.na.t), "missing values in the time column\n")
        if (sum(lst.na.x) > 0) 
            cat("\n\nThere are ", sum(lst.na.x), "missing values in the x column(s)\n")
        if (sum(lst.na.y) > 0) 
            cat("\n\nThere are ", sum(lst.na.y), "missing values in the y column\n")
        if (sum(lst.na.sd) > 0) 
            cat("\n\nThere are ", sum(lst.na.sd), "missing values in the sd column\n")
        if (sum(lst.na.nn) > 0) 
            cat("\n\nThere are ", sum(lst.na.nn), "missing values in the group size column\n")
        if (sum(lst.na.covar) > 0) 
            cat("\n\nThere are ", sum(lst.na.covar), "missing values in the covariate column\n")
        if (sum(lst.na.fct1) > 0) 
            cat("\n\nThere are ", sum(lst.na.fct1), "missing values in the covariate column on a\n")
        if (sum(lst.na.fct2) > 0) 
            cat("\n\nThere are ", sum(lst.na.fct2), "missing values in the covariate column in b\n")
        if (sum(lst.na.fct3) > 0) 
            cat("\n\nThere are ", sum(lst.na.fct3), "missing values in the covariate column on var\n")
        if (sum(lst.na.fct4) > 0) 
            cat("\n\nThere are ", sum(lst.na.fct4), "missing values in the covariate column on c\n")
        if (sum(lst.na.fct5) > 0) 
            cat("\n\nThere are ", sum(lst.na.fct5), "missing values in the covariate column on d\n")
    }
    lst.na <- (lst.na.x | lst.na.t | lst.na.y | lst.na.sd | lst.na.nn | 
        lst.na.covar | lst.na.fct1 | lst.na.fct2 | lst.na.fct3 | 
        lst.na.fct4 | lst.na.fct5)
    if (output) 
        if (sum(lst.na) > 0) 
            cat("\n these rows are removed from the dataset\n")
    dfr <- dfr[!lst.na, ]
    if (exists("track")) 
        print("f.remove.NAs:  END")
    return(dfr)
}



f.cont <- function(ans.all) {
    if (exists("track")) 
        print("f.cont")
    if (ans.all$dtype %in% c(1, 5, 15, 10, 25, 26, 250, 260)) 
        cont <- T
    if (ans.all$dtype %in% c(2, 3, 4, 6, 84)) 
        cont <- F
    return(cont)
    print("f.cont: END")
}



f.check.nonneg.num <- function(vec, gui = FALSE, dtype = NA, quick.ans = 1) {
    if (exists("track")) 
        print("f.check.nonneg.num")
    check.out <- F
    if (quick.ans == 1) {
        if (!is.numeric(vec)) {
            cat("\na problem was encountered in a column of the data\n")
            cat("the values are:")
            print(vec)
            cat("\nif you do not see any anomaly, check if there might be blanks in this column\n")
            cat("e.g. by using the filter option in Excel\n\n")
            stop("\nATTENTION: not all values are numerical\n")
        }
        if (any(vec < 0)) {
            if (gui) 
                stop("\nvariable contain negative values\n", 
                  call. = F)
            if (!gui) {
                cat("\n\nATTENTION: variable contains negative values\n")
                if (!(dtype %in% c(25, 250))) {
                  ans.cont <- menu(c("no, this is an error, stop analysis", 
                    "yes, but omit log-transformation"), title = "\nDo you want to continue nonetheless?")
                  switch(ans.cont, stop(message = "PROAST interrupted", 
                    call. = F), check.out <- T)
                }
            }
        }
        if (!any(vec > 0)) {
            if (gui) 
                stop("\nvariable contains zeros only\n", call. = F)
            if (!gui) {
                cat("\n\nATTENTION: variable contains zeros only\n")
                ans.cont <- 0
                switch(ans.cont, stop(message = "PROAST interrupted", 
                  call. = F), check.out <- T)
            }
        }
    }
    if (quick.ans > 1) {
        if (!is.numeric(vec) || any(vec < 0) || !any(vec > 0)) {
            cat("\nATTENTION: response data contain negative values, or only zeros,\n              or non-numerical values\nno analysis possible")
            print("the response data are:")
            print(vec)
            check.out <- T
        }
    }
    if (exists("track")) 
        print("f.check.nonneg.num:   END")
    invisible(check.out)
}



f.execute <- function(ans.all, no.plot = F) {
    if (exists("track")) 
        print("f.execute")
    ans.all$no.plot <- no.plot
    ans.all$cont <- f.cont(ans.all)
    with(ans.all, {
        if (length(xans) > 1) {
            if (ans.all$aic.crit == 2) 
                ans.all$aic.crit <- 5
            if (model.ans != 42) 
                x.leg <- "dose after addition"
            if (model.ans != 42) 
                if (sf.x != 1) 
                  x.leg <- paste("dose/", sf.x, " after addition", 
                    sep = "")
            Mx <- data.0[, xans]/sf.x
            ans.all$full.ans <- 2
        }
        if (length(xans) == 1) {
            x <- as.numeric(data.0[, xans])/sf.x
            f.check.nonneg.num(x, gui)
            x.leg <- varnames[xans]
            if (sf.x != 1) 
                x.leg <- paste(x.leg, "/", sf.x, sep = "")
        }
        if (length(xans) > 1) {
            data.0 <- data.0[order(Mx[, 1]), ]
            Mx <- data.0[, xans]/sf.x
            x <- rowSums(Mx)/sf.x
        }
        else {
            data.0 <- data.0[order(x), ]
            x <- as.numeric(data.0[, xans])/sf.x
        }
        if (tans > 0) {
            ttt <- as.numeric(data.0[, tans])/sf.t
            f.check.nonneg.num(ttt, gui)
            t.leg <- varnames[tans]
            ttt.txt <- levels(factor(ttt))
        }
        if (cont) {
            if (yans != 0) {
                y <- as.numeric(data.0[, yans])/sf.y
                yy <- y
                y.leg <- varnames[yans]
                if (sf.y != 1) 
                  y.leg <- paste(y.leg, "/", sf.y, sep = "")
                if (dtype %in% c(1, 5, 10, 15)) {
                  data.problem <- f.check.nonneg.num(y, gui, 
                    dtype)
                  if (data.problem) {
                    cat("\nATTENTION: analyses will be done without log-transformation, due to negative observations\n\n")
                    if (dtype == 1) 
                      dtype <- 25
                    if (dtype %in% c(10)) 
                      dtype <- 250
                  }
                }
            }
            if (auto.detlim) {
                if (min(y) == 0) 
                  detlim <- min(y[y > 0])
                detlim.col <- nvar + 1
                data.0[, detlim.col] <- detlim
                Vdetlim <- data.0[, detlim.col]
            }
            else if (detlim.col > 0) 
                Vdetlim <- data.0[, detlim.col]
            else if (detlim > 0) {
                Vdetlim <- rep(detlim, length(data.0[, 1]))
                data.0 <- cbind(data.0, Vdetlim)
            }
            if (dtype %in% c(1, 5) && min(y) == 0 && min(Vdetlim) == 
                0) {
                cat("\n\nATTENTION:  minimum value of observed response data is zero\n")
                if (dtype %in% c(10, 15, 250, 260)) 
                  cat("(note that detection limit does not really apply to summary data)\n")
                ans.detlim <- menu(c("Give detection limit", 
                  "Add nonzero value to observations"), title = "Q6-7: What do you want to do?")
                if (ans.detlim == 1) {
                  cat("\nthe lowest nonzero value is: ", min(y[y > 
                    0]), "\n")
                  detlim <- eval(parse(prompt = paste("\n\nQ6: Give global value of detection limit,\n                       or type 0 if you have individual values > ")))
                  if (detlim == 0) {
                    detlim.col <- menu(varnames[1:nvar], title = "\nQ6a: Give column with individual detection limits\n")
                    Vdetlim <- data.0[, detlim.col]
                  }
                  if (detlim != 0 && detlim.col > 0) {
                    Vdetlim <- data.0[, detlim.col]
                    if (min(Vdetlim <= 0)) {
                      cat("\n ATTENTION: there are non-positive detection limits !\n")
                      f.press.key.to.continue()
                    }
                    if (any(is.na(Vdetlim))) {
                      cat("\nATTENTION:  column with detection limits contains NAs; \n                       replace NAs with (large) numbers\n\n")
                      f.press.key.to.continue()
                    }
                  }
                  if (detlim.col == 0) {
                    Vdetlim <- rep(detlim, length(y))
                    data.0 <- cbind(data.0, Vdetlim)
                    detlim.col <- nvar + 1
                  }
                }
                if (ans.detlim == 2) {
                  cat("\nThe lowest nonzero observation is: ", 
                    min(y[y > 0]), "\n")
                  y.add <- eval(parse(prompt = "\n\nQ7: Give small value to be added to the data > "))
                }
            }
            if (min(Vdetlim) > 0) {
                y[y < Vdetlim] <- -1000
                yy <- yy * (y != -1000) + (0.5 * Vdetlim) * (y == 
                  -1000)
            }
            if (!is.na(cens.up)) {
                y[y >= cens.up] <- -2000
                yy <- yy * (y != -2000) + (1 * cens.up) * (y == 
                  -2000)
            }
            if (max(Vdetlim) > 0) 
                ans.all$low.y <- min(Vdetlim)/2
            else ans.all$low.y <- 0.98 * min(yy)
            ans.all$upp.y <- 1.02 * max(yy)
            if (ans.detlim == 2) {
                y <- y + y.add
                yy <- y
            }
            if (dtype %in% c(1, 5)) {
                y.denom <- 1
                if (yans.denom != 0) {
                  y.denom <- data.0[, yans.denom]
                  if (min(y.denom) == 0) 
                    cat("\nATTENTION: zero values in denominator ! \n\n")
                  y <- y/y.denom
                  yy <- yy/y.denom
                  y.leg <- paste(y.leg, "/", varnames[yans.denom])
                }
            }
            nn <- 0
            sd <- NA
            if (dtype %in% c(10, 15, 250, 260)) {
                if (sans == 0) 
                  stop(message = "no data for variation statistics provided\n")
                if (nans == 0) 
                  stop(message = "no data for group sizes provided\n")
                sd <- as.numeric(data.0[, sans])/sf.y
                nn <- as.numeric(data.0[, nans])
                if (any(y == 0)) {
                  cat("\nATTENTION: there are zero (mean) observations in the data\n")
                  cat("\nThis is not allowed for mean values\n")
                  cat("the analysis will be interrupted\n")
                  stop()
                }
                if (any(sd == 0)) {
                  cat("\nATTENTION: there are zero SDs in the data\n")
                  cat("\nThis is not allowed for SDs\n")
                  cat("the analysis will be interrupted\n")
                  stop()
                }
                if (sd.se == 2) 
                  sd <- sd * sqrt(nn)
                if (0) {
                  f.graph.window(2)
                  plot(log10(y), log10(sd), xlab = "mean (log10-scale)", 
                    ylab = "sd (log10-scale)")
                  title(main = "correlation between mean of y and sd of y")
                  y.mn <- y
                  CV <- sd/y.mn
                  mn.log <- logb(y.mn/sqrt(1 + CV^2))
                  mn.log <- mn.log/logb(10)
                  sd2.log <- logb(CV^2 + 1)
                  sd2.log <- sd2.log/logb(10)^2
                  plot(mn.log, sd2.log, xlab = "mean of log-y", 
                    ylab = "sd of log-y")
                  title(main = "correlation between mean of log-y and sd of log-y")
                  eval(parse(prompt = "\n\nCheck plot and type c to continue .................. >   "))
                }
                if (dtype %in% c(10, 15)) {
                  y.mn <- y
                  CV <- sd/y.mn
                  mn.log <- logb(y.mn/sqrt(1 + CV^2))
                  yy <- exp(mn.log)
                  ans.all$mn.log <- mn.log
                  ans.all$sd2.log <- logb(CV^2 + 1)
                }
                if (dtype == 250) {
                  y.mn <- y
                  ans.all$mn.log <- y.mn
                  ans.all$sd2.log <- sd^2
                  ans.all$yy <- y.mn
                }
            }
        }
        if (dtype %in% c(4, 6, 84)) {
            nn <- data.0[, nans]
            if (yans != 0) {
                kk <- data.0[, yans]
                nn <- data.0[, nans]
                y.leg <- varnames[yans]
            }
            if (sum(kk > nn) > 0) {
                cat("\nNumber of responses is larger than sample size ! \n\n")
                stop()
            }
            y <- kk/nn
            yy <- y
            if (max(y, na.rm = T) < 0.5) 
                ans.all$y.lim.plt <- c(0, 1.2 * max(y))
            else ans.all$y.lim.plt <- c(0, 1)
        }
        if (!cont) {
            cens <- 0
            if (sum(varnames == "cens") > 0) {
                zzz <- 1:nvar
                cens.col <- zzz[varnames == "cens"]
                cens <- data.0[, cens.col]
            }
        }
        if (dtype %in% 2:3) {
            if (yans != 0) {
                y <- as.numeric(data.0[, yans])
                yy <- y
                y.leg <- varnames[yans]
            }
            nn <- 1
            y.original <- y
            scores.orig <- sort(unique(y))
            ymn <- tapply(y, x, mean)
            xmn <- tapply(x, x, mean)
            ymn <- as.numeric(ymn)
            if (reverse.ans == 0) {
                if (length(xmn) > 1) 
                  if (var(xmn != 0) && var(ymn != 0)) 
                    if (cor(xmn, ymn) < 0) {
                      cat("\nProast assumes zero to represent normal, and higher scores abnormal")
                      cat("\nYour data seem to have the opposite direction\n")
                      ans.all$reverse.ans <- menu(c("yes", "no"), 
                        title = "Do your scores have opposite direction?\n")
                    }
            }
            if (ans.all$reverse.ans == 1) 
                scores.orig <- rev(scores.orig)
            score <- 0
            n.lev <- length(scores.orig)
            for (ii in 1:n.lev) {
                y[y.original == scores.orig[ii]] <- score
                score <- score + 1
            }
            scores.mtr <- cbind(scores.orig, levels(factor(y)))
            dimnames(scores.mtr) <- list(NULL, c("orig.scores", 
                "temp.scores"))
            dum.ord <- sum(scores.mtr[, 1] == scores.mtr[, 2]) != 
                length(scores.mtr[, 1])
            if (dum.ord) {
                cat("\n\nATTENTION: the original scores have been transformed for analysis as follows:\n")
                print(scores.mtr)
            }
            ans.all$y.transf <- y
        }
        fct1.txt <- ""
        fct2.txt <- ""
        fct3.txt <- ""
        fct4.txt <- ""
        fct5.txt <- ""
        fct1 <- rep(1, length(y))
        fct2 <- rep(1, length(y))
        fct3 <- rep(1, length(y))
        fct4 <- rep(1, length(y))
        fct5 <- rep(1, length(y))
        if (quick.ans == 1) {
            if (fct1.no > 0) {
                fct1 <- as.numeric(factor(data.0[, fct1.no]))
                fct1.txt <- levels(factor(data.0[, fct1.no]))
                if (interrupt && max(fct1) == 1) {
                  cat("\nAttention 1: The factor you chose as covariate on par a has only one level\n")
                  cat("you might have selected a subgroup for this factor\n")
                  f.press.key.to.continue()
                }
            }
            if (fct2.no > 0) {
                fct2 <- as.numeric(factor(data.0[, fct2.no]))
                fct2.txt <- levels(factor(data.0[, fct2.no]))
                if ((cont || model.type == 2) && model.ans %in% 
                  c(32, 45)) {
                  cat("\nATTENTION: covariate on parameter b is not possible in these NMDRs\n\n")
                  fct2 <- 1
                }
                if (max(fct2) == 1) {
                  cat("\nAttention: the factor you chose as covariate on par b has only one level\n")
                  cat("you might have selected a subgroup for this factor\n")
                  f.press.key.to.continue()
                }
            }
            if (fct3.no > 0) {
                fct3 <- data.0[, fct3.no]
                fct3.txt <- levels(factor(data.0[, fct3.no]))
                if (cont) {
                  fct3 <- as.numeric(factor(data.0[, fct3.no]))
                  if (max(fct3) == 1) {
                    cat("\nAttention: the factor you chose as covariate on var has only one level\n")
                    cat("you might have selected a subgroup for this factor\n")
                    f.press.key.to.continue()
                  }
                }
                else if (dtype == 4) {
                  fct3 <- log(data.0[, fct3.no])
                  ans.all$th.par <- log(as.numeric(levels(as.factor(data.0[, 
                    fct3.no]))))
                }
                else if (dtype == 6) {
                  fct3 <- as.numeric(factor(data.0[, fct3.no]))
                }
                ans.all$nr.var <- max(fct3)
            }
            else ans.all$nr.var <- 1
            if (fct4.no > 0) {
                fct4 <- as.numeric(factor(data.0[, fct4.no]))
                fct4.txt <- levels(factor(data.0[, fct4.no]))
                if ((cont || model.type == 2) && model.ans %in% 
                  c(32, 45)) {
                  cat("\nATTENTION: covariate on parameter c is not possible in NMDRs\n\n")
                  fct4 <- 1
                }
                if (max(fct4) == 1) {
                  cat("\nAttention: the factor you chose as covariate on par c has only one level\n")
                  cat("you might have selected a subgroup for this factor\n")
                  f.press.key.to.continue()
                }
            }
            if (fct5.no > 0) {
                fct5 <- as.numeric(factor(data.0[, fct5.no]))
                fct5.txt <- levels(factor(data.0[, fct5.no]))
                if ((cont || model.type == 2) && model.ans %in% 
                  c(31, 32, 33, 45)) {
                  cat("\nATTENTION: covariate on parameter d is not possible in NMDRs\n\n")
                  fct5 <- 1
                }
                if (max(fct5) == 1) {
                  cat("\nAttention: the factor you chose as covariate on par d has only one level\n")
                  cat("you might have selected a subgroup for this factor\n")
                  f.press.key.to.continue()
                }
            }
            twice <- F
            if (fct1.no != 0) 
                twice <- fct1.no == fct2.no
            if (!twice && fct1.no > 1) 
                twice <- fct1.no == fct5.no
            if (!twice && fct1.no > 1) 
                twice <- fct1.no == fct4.no
        }
        if (quick.ans > 1) {
            twice <- T
            covariate <- rep(1, length(y))
            covar.txt <- ""
            if (covar.no > 0) {
                covariate <- as.numeric(factor(data.0[, covar.no]))
                if (max(covariate, na.rm = T) == 1) {
                  cat("\nAttention: the factor you chose as covariate has only one level\n")
                  ans.all$covar.no <- 0
                }
                covar.txt <- levels(factor(data.0[, covar.no]))
                covar.name <- varnames[covar.no]
                if (sum(is.na(covariate)) > 0) {
                  cat("\nAttention: the factor you chose as covariate contains NAs\n")
                  cat("\nThe associated observations will be removed\n")
                  f.press.key.to.continue()
                  covariate <- data.0[, covar.no]
                }
                fct1.txt <- covar.txt
                fct2.txt <- covar.txt
                if (cont) {
                  fct3.txt <- covar.txt
                }
            }
            ans.all$covariate <- covariate
            ans.all$covar.txt <- covar.txt
            ans.all$covar.name <- covar.name
        }
        ans.all$nr.aa <- max(fct1)
        ans.all$nr.bb <- max(fct2)
        ans.all$nr.cc <- max(fct4)
        ans.all$nr.dd <- max(fct5)
        ans.all$fct1 <- fct1
        ans.all$fct2 <- fct2
        ans.all$fct3 <- fct3
        ans.all$fct4 <- fct4
        ans.all$fct5 <- fct5
        ans.all$fct1.txt <- fct1.txt
        ans.all$fct2.txt <- fct2.txt
        if (length(xans) > 1) 
            ans.all$fct2.txt <- substring(varnames[xans[1:length(xans)]], 
                1, 3)
        ans.all$fct3.txt <- fct3.txt
        ans.all$fct4.txt <- fct4.txt
        ans.all$fct5.txt <- fct5.txt
        ans.all$cont <- cont
        if (0) 
            if (length(xans) == 1) {
                if (quick.ans == 1) {
                  if (length(unique(fct3)) > 1 && dtype == 4) {
                    ans.all$get.gr.txt <- TRUE
                    ans.all$gr.txt <- f.pars.frq(ans.all)$gr.txt
                  }
                  else ans.all$gr.txt <- f.pars(ans.all)$gr.txt
                }
                if (quick.ans > 1) 
                  ans.all$gr.txt <- covar.txt
            }
        if (displ.no > 0) {
            ans.all$displ.fact <- factor(data.0[, displ.no])
            ans.all$displ.txt <- levels(ans.all$displ.fact)
            if (model.ans == 42) 
                ans.all$fct2 <- ans.all$displ.fact
        }
        if (displ.no == 0) {
            ans.all$displ.fact <- 1
            ans.all$displ.txt <- NA
        }
        ans.all$displ.no <- displ.no
        if (quick.ans == 1) 
            trace.plt <- T
        else trace.plt <- F
        if (quick.ans == 1) 
            trace <- T
        else trace <- F
        if (min(x) == 0) 
            dum.contr <- min(x[x > 0])/4
        else dum.contr <- min(x)
        if (min(y) == 0) 
            dum.zero.resp <- min(yy[yy > 0])/2
        else dum.zero.resp <- min(yy)
        xy.lim[1] <- dum.contr
        xy.lim[2] <- min(x)
        if (xy.lim[2] > 0) {
            xy.lim[2] <- 0
            xy.lim[1] <- min(x)/3
        }
        xy.lim[3] <- max(x)
        if (0) 
            if (length(xans) > 1) {
                xy.lim[2] <- min(Mx)
                xy.lim[3] <- max(Mx)
            }
        xy.lim[4] <- min(yy)
        xy.lim[5] <- max(yy)
        if (!cont) {
            xy.lim[4] <- 0
            xy.lim[5] <- 1
            if (max(y) < 0.5) 
                xy.lim[5] <- (max(y) + 0.1) * 1.2
        }
        xy.lim[6] <- dum.zero.resp
        xy.lim[3] <- 1.05 * xy.lim[3]
        xy.lim[5] <- 1.05 * xy.lim[5]
        xy.lim <- signif(xy.lim, 3)
        if (WAPP) 
            ans.all$logprob <- FALSE
        if (cont) {
            if (plt.mns == 2) 
                cex.1 <- 1
            else cex.1 <- 0.75
            if (plot.type == 9) 
                cex.1 <- 0.6
            if (dtype %in% c(10, 15, 250, 260)) 
                cex.1 <- 1.6
            cex.2 <- 1.2
            if (dtype == 5) 
                cex.1 <- 0.5
            if (dtype == 15) 
                cex.2 <- 0.75
        }
        if (!cont) {
            if (dtype == 6) 
                cex.1 <- 0.6
            else cex.1 <- 1.2
            cex.2 <- 1.5
        }
        ans.all$cex.1 <- cex.1
        ans.all$cex.2 <- cex.2
        ans.all$quick.ans <- quick.ans
        if (quick.ans > 1) 
            ans.all$quick <- T
        ans.all$data.0 <- data.0
        ans.all$x <- x
        if (length(xans) > 1) {
            ans.all$nr.dosecol <- length(Mx[1, ])
            ans.all$Mx <- as.matrix(Mx)
        }
        else {
            ans.all$nr.dosecol <- 1
            ans.all$Mx <- NA
        }
        ans.all$y <- y
        ans.all$yy <- yy
        ans.all$dtype <- dtype
        ans.all$cont <- cont
        if (cont) 
            ans.all$model.names <- f.expect.con(name = T)
        ans.all$x.leg <- x.leg
        if (yans > 0) 
            ans.all$y.leg <- y.leg
        ans.all$sd <- sd
        ans.all$nn <- nn
        ans.all$detlim <- detlim
        ans.all$Vdetlim <- Vdetlim
        ans.all$detlim.col <- detlim.col
        ans.all$ans.detlim <- ans.detlim
        ans.all$y.add <- y.add
        if (dtype %in% c(4, 6, 84) && yans != 0) {
            ans.all$kk <- kk
            ans.all$y <- kk/nn
            ans.all$cens <- cens
        }
        ans.all$ttt <- ttt
        ans.all$ttt.txt <- ttt.txt
        ans.all$t.leg <- t.leg
        t
        ans.all$xy.lim <- xy.lim
        ans.all$twice <- twice
        ans.all$factor.name <- ""
        ans.all$factor.name <- f.factorname(ans.all)
        ans.all$low.y <- low.y
        ans.all$upp.y <- upp.y
        if (dtype == 11) 
            ans.all$fit.ans <- 2
        if (ans.all$fit.ans == 2) {
            if (!cont) 
                print("fit.ans has value 2 !!")
            ans.all$plot.type <- 4
        }
        ans.all$trace <- trace
        ans.all$trace.plt <- trace.plt
        ans.all$created <- date()
        if (dtype == 3) {
            ans.all$model.type <- 2
            ans.all$ces.ans <- 1
            ans.all$scores.orig <- scores.orig
            ans.all$scores.mtr <- scores.mtr
            ans.all$nth <- max(y)
            ans.all$plot.type <- 6
            if (ans.all$cat.ans > 0) {
                y[y < cat.ans] <- 0
                y[y >= cat.ans] <- 1
                ans.all$dtype <- 2
                ans.all$y <- y
                ans.all$nth <- 1
                ans.all$plot.type <- 2
            }
        }
        if (dtype %in% c(2, 3, 6, 15)) 
            ans.all$CI.plt <- FALSE
        if (dtype %in% c(4, 10, 250)) 
            ans.all$CI.plt <- TRUE
        if (!no.plot) {
            f.graph.window(1)
            if (!fitted) 
                ans.all$heading <- "data"
            if (cont) 
                f.plot.con(ans.all)
            if (dtype %in% c(4, 6, 84)) {
                if (max(fct2) == 1) 
                  kk.dum <- kk
                else kk.dum <- 1
                if (dtype == 6) 
                  kk.dum <- 1
                if (length(y) > 200) 
                  kk.dum <- 1
                if (ans.all$plot.type == 0) 
                  ans.all$plot.type <- 1
                f.plot.frq(ans.all)
            }
            if (dtype %in% 2:3) 
                f.plot.cat(plot.type = 1, x, y, x.leg, y.leg, 
                  fct1 = fct1, fct2 = fct2, shift = 0, heading = "", 
                  color = color)
        }
        if (dtype != 6 && !do.MA) {
            ans.all$pi.full <- NA
            ans.all$x.full <- NA
            ans.all$fct1.full <- NA
            ans.all$fct2.full <- NA
            ans.all$loglik.full <- NA
            ans.all$npar.full <- NA
        }
        if (dtype == 6) {
            ans.all$nr.alfa <- max(fct3)
        }
        if (0) 
            if (!gui) {
                dtype.tmp <- odt$dtype[ans.all$yans]
                if (dtype.tmp != ans.all$dtype) {
                  cat("\nThe type of data you indicated does not match the one indicated in the dataset\n")
                  cat("indicated in dataset:", dtype.tmp)
                  cat("\nintended:", ans.all$dtype, "\n")
                  f.press.key.to.continue()
                }
            }
        if (cont && length(xans) == 1) 
            ans.all$x.mn <- mean(x)
        if (cont) {
            dum.nn <- f.nr.replicates(ans.all)
            if (mean(dum.nn) < 3) 
                ans.all$plt.mns <- 2
            if (length(xans) > 1) {
                ans.all$plt.mns <- 3
                ans.all$CI.plt <- TRUE
            }
        }
        ans.all$dtype.0 <- dtype
        if (lump.ans == 1) 
            ans.all <- f.lump.bin(ans.all)
        if (incr.decr.no > 0) {
            ans.all$fct4 <- as.numeric(factor(data.0[, incr.decr.no]))
            ans.all$nr.cc <- 2
            ans.all$fct4.txt <- c("decr", "incr")
        }
        if (quick.ans %in% c(1, 6)) 
            ans.all$model.fam <- 1
        else ans.all$model.fam <- 0
        f.assign(".Pr.last", ans.all)
        if (exists("track")) 
            print("f.execute:  END")
        return(ans.all)
    })
}



f.clear <- function(ans.all) {
    if (exists("track")) 
        print("f.clear")
    ans.all$show <- ""
    ans.all$MLE <- NA
    ans.all$regr.par <- NA
    ans.all$regr.par.matr <- NA
    ans.all$CED.matr <- matrix(NA)
    ans.all$regr.par.matr <- matrix(NA)
    ans.all$conf.int <- matrix(NA, ncol = 2)
    ans.all$fitted <- F
    ans.all$ans.plt <- 3
    ans.all$par.start <- NA
    ans.all$scale.dum <- NA
    ans.all$lb <- NA
    ans.all$ub <- NA
    ans.all$lower <- NA
    ans.all$upper <- NA
    ans.all$loglik <- NA
    ans.all$do.MA <- FALSE
    ans.all$MA.running <- FALSE
    if (ans.all$quick.ans > 1) {
        ans.all$Vbmdl <- NA
        ans.all$Vbmdu <- NA
    }
    if (exists("track")) 
        print("f.clear:   END")
    return(ans.all)
}



f.full.ans <- function(ans.all, gui) {
    if (exists("track")) 
        print("f.full.ans")
    if (gui == T) 
        ans.all$gui <- T
    with(ans.all, {
        if (dtype %in% c(2, 3)) {
            full.ans <- 2
            return(full.ans)
        }
        if (dtype %in% c(6)) {
            full.ans <- 1
            return(full.ans)
        }
        dum.nn <- f.nr.replicates(ans.all)
        if (cont) {
            ans.all$model.ans <- 11
            regr.par <- f.regr.par.full(ans.all)
            nr.dosegr <- length(regr.par)
        }
        if (!cont) {
            ans.all$model.type <- 1
            ans.all$model.ans <- 14
            ans.all <- f.start.bin(ans.all, tmp.quick = T)
            nr.dosegr <- length(ans.all$regr.par)
        }
        full.ans <- ans.all$full.ans
        if (!gui && dtype != 15) {
            if (mean(dum.nn) < 3) {
                cat("\n\naverage group size is only", round(mean(dum.nn), 
                  1))
                full.ans <- menu(c("yes", "no"), title = "\nDo you nonetheless want to fit full model?\n")
            }
            else if (nr.dosegr > 100 || length(dum.nn) > 100) {
                cat("\ndataset is large (total number of dose groups: ", 
                  nr.dosegr, "), fitting full model may take long")
                full.ans <- menu(c("yes", "no"), title = "\nDo you nonetheless want to fit full model?\n")
            }
        }
        if (gui) {
            full.ans <- 1
            if (mean(dum.nn) < 3) {
                cat("\n\naverage group size is only", round(mean(dum.nn), 
                  1))
                cat("\ntherefore full model will not be fitted\n")
                full.ans <- 2
            }
            else if (nr.dosegr > 100 || length(dum.nn) > 100) {
                cat("\ndataset is large, therefore full model will not be fitted\n")
                full.ans <- 2
            }
        }
        if (exists("track")) 
            print("f.full.ans:   END")
        return(full.ans)
    })
}



f.regr.par.full <- function(ans.all) {
    if (exists("track")) 
        print("f.regr.par.full")
    with(ans.all, {
        if (dtype %in% c(10, 15, 250, 260)) {
            if (dtype %in% c(10, 15)) 
                regr.par <- exp(mn.log)
            if (dtype == 250) 
                regr.par <- (mn.log)
            if (dtype == 260) 
                regr.par <- (mn.log)^2
        }
        else {
            regr.par <- numeric()
            nn <- numeric()
            if (length(xans) > 1) 
                x <- data.0[, xans[1]]
            if (nr.var == 1) {
                if (twice) 
                  for (jj in levels(factor(fct2))) {
                    x.tmp <- x[fct2 == jj]
                    y.tmp <- yy[fct2 == jj]
                    regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                      x.tmp, mean)))
                    regr.par <- c(regr.par, regr.tmp)
                    nn.tmp <- tapply(y.tmp, x.tmp, length)
                    nn <- c(nn, nn.tmp)
                  }
                if (!twice) 
                  for (jj in levels(factor(fct2))) for (ii in levels(factor(fct1))) {
                    x.tmp <- x[fct1 == ii & fct2 == jj]
                    y.tmp <- yy[fct1 == ii & fct2 == jj]
                    regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                      x.tmp, mean)))
                    regr.par <- c(regr.par, regr.tmp)
                    nn.tmp <- tapply(y.tmp, x.tmp, length)
                    nn <- c(nn, nn.tmp)
                  }
            }
            if (nr.var > 1) {
                regr.par <- numeric()
                for (jj in levels(factor(fct3))) {
                  x.tmp <- x[fct3 == jj]
                  y.tmp <- yy[fct3 == jj]
                  regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                    x.tmp, mean)))
                  regr.par <- c(regr.par, regr.tmp)
                  nn.tmp <- tapply(y.tmp, x.tmp, length)
                  nn <- c(nn, nn.tmp)
                }
            }
            ans.all$nn <- nn
        }
        return(regr.par)
    })
}



f.constr.dd <- function(model.ans) {
    if (exists("track")) 
        print("f.constr.dd")
    lower.dd <- 0.2
    upper.dd <- 5
    if (model.ans %in% c(18, 20, 21, 23, 25)) {
        upper.dd <- upper.dd * 1.5
    }
    if (model.ans %in% 51:52) {
        lower.dd <- lower.dd/5
    }
    if (model.ans %in% 53:54) {
        lower.dd <- lower.dd/4
    }
    if (exists("track")) 
        print("f.constr.dd:    END")
    return(c(lower.dd, upper.dd))
}



f.select.con <- function(ans.all, interactive_mode = T, results_env = NULL, display_plots = TRUE) {
    if (exists("track")) 
        print("f.select.con")
    ans.all$HILL <- NA
    ans.all$INVEXP <- NA
    ans.all$LOGN <- NA
    if (ans.all$nr.models == 1) 
        ans.all$model.fam <- 1
    else ans.all$model.fam <- 0
    if (ans.all$nr.models == 2) 
        ans.all$model.list <- c("EXP", "HILL")
    if (ans.all$nr.models == 3) 
        ans.all$model.list <- c("EXP", "HILL", "INVEXP")
    if (ans.all$nr.models == 4) 
        ans.all$model.list <- c("EXP", "HILL", "INVEXP", "LOGN")
    Vaic <- numeric(ans.all$nr.models)
    Vnpar.sel <- numeric()
    if (ans.all$output) 
        cat(paste("\n\nresponse: ", ans.all$y.leg))
    ans.all$model.switch <- 0
    if (ans.all$quick.ans == 6) {
        ans.all <- f.select.m46.con(ans.all, interactive_mode = interactive_mode, results_env = results_env, display_plots = display_plots)
        return(ans.all)
        # assign(ans.all$created, ans.all, envir = results_env)
    }
    if (ans.all$NES.ans == 2) 
        ans.all$CES <- 0.05
    ans.all <- f.select.m5.con(ans.all, output = ans.all$output, interactive_mode = interactive_mode, results_env = results_env, display_plots = display_plots)
    Vaic[1] <- ans.all$EXP$aic
    Vnpar.sel <- c(Vnpar.sel, ans.all$EXP$npar)
    ans.all$TREND <- TRUE
    if (!ans.all$EXP$trend) 
        ans.all$TREND <- FALSE
    if (ans.all$nr.models > 1) {
        if (ans.all$output) 
            cat(paste("\n\nresponse: ", ans.all$y.leg))
        ans.all$model.switch <- 1
        ans.all <- f.select.m5.con(ans.all, output = ans.all$output, interactive_mode = interactive_mode, results_env = results_env, display_plots = display_plots)
        Vaic[2] <- ans.all$HILL$aic
        Vnpar.sel <- c(Vnpar.sel, ans.all$HILL$npar)
        if (!ans.all$TREND) 
            if (ans.all$HILL$trend) 
                ans.all$TREND <- TRUE
    }
    if (ans.all$nr.models > 2) {
        if (ans.all$output) 
            cat(paste("\n\nresponse: ", ans.all$y.leg))
        ans.all$model.switch <- 2
        ans.all <- f.select.m5.con(ans.all, output = ans.all$output, interactive_mode = interactive_mode, results_env = results_env, display_plots = display_plots)
        Vaic[3] <- ans.all$INVEXP$aic
        Vnpar.sel <- c(Vnpar.sel, ans.all$INVEXP$npar)
        if (!ans.all$TREND) 
            if (ans.all$INVEXP$trend) 
                ans.all$TREND <- TRUE
    }
    if (ans.all$nr.models > 3) {
        if (ans.all$output) 
            cat(paste("\n\nresponse: ", ans.all$y.leg))
        ans.all$model.switch <- 3
        ans.all <- f.select.m5.con(ans.all, output = ans.all$output, interactive_mode = interactive_mode, results_env = results_env, display_plots = display_plots)
        Vaic[4] <- ans.all$LOGN$aic
        Vnpar.sel <- c(Vnpar.sel, ans.all$LOGN$npar)
        if (!ans.all$TREND) 
            if (ans.all$LOGN$trend) 
                ans.all$TREND <- TRUE
    }
    aic.min <- min(Vaic)
    model.sel <- (1:ans.all$nr.models)[Vaic == aic.min]
    if (0) 
        if (ans.all$full.ans == 1) {
            cat("\nAIC of full model = ", ans.all$aic.full, "\n")
            cat("AIC of best fitting model = ", aic.min, "\n")
        }
    npar.sel <- max(Vnpar.sel[model.sel])
    if (length(ans.all$xans) == 1 && !ans.all$model.ans %in% 
        c(5, 15)) 
        if (ans.all$full.ans == 1) 
            if (ans.all$npar.aic.full >= npar.sel) 
                if (aic.min > ans.all$aic.full + ans.all$aic.crit) {
                  cat("\nAIC of full model = ", ans.all$aic.full, 
                    "\n")
                  cat("AIC of best fitting model = ", aic.min, 
                    "\n")
                  ans.all$alert.full <- f.alert.full()
                  if (ans.all$output) 
                    cat(ans.all$alert.full)
                }
    ans.all$Vaic <- Vaic
    ans.all$fitted <- T
    if (exists("track")) 
        print("f.select.con END")
    return(ans.all)
    #assign(ans.all$created, ans.all, envir = results_env)
}



f.select.m5.con <- function(ans.all, output = TRUE, interactive_mode = T, results_env = NULL, display_plots = TRUE) {
    if (exists("track")) 
        print("f.select.m5.con")
    gui <- ans.all$gui
    model.switch <- ans.all$model.switch
    cont <- ans.all$cont
    quick.ans <- ans.all$quick.ans
    x <- ans.all$x
    y <- ans.all$y
    yy <- ans.all$yy
    nn <- ans.all$nn
    if (ans.all$dtype.0 == 5) 
        ans.all$dtype <- 5
    if (ans.all$dtype.0 == 15) 
        ans.all$dtype <- 15
    dtype <- ans.all$dtype
    if (dtype == 5) 
        ans.all$dtype.0 <- 5
    nth <- ans.all$nth
    covariate <- ans.all$covariate
    covar.no <- ans.all$covar.no
    covar.txt <- ans.all$covar.txt
    CES <- ans.all$CES
    ces.ans <- ans.all$ces.ans
    CED <- NA
    xans <- ans.all$xans
    cens.up <- ans.all$cens.up
    select.name <- ans.all$select.name
    factor.name <- ans.all$factor.name
    if (cont) 
        plot.type <- ans.all$plot.type
    else plot.type <- 0
    color <- ans.all$color
    if (dtype == 6) 
        alfa.length <- 1
    else alfa.length <- 0
    alfa.start <- ans.all$alfa.start
    full.ans <- ans.all$full.ans
    interrupt <- ans.all$interrupt
    nr.covar <- max(covariate, na.rm = TRUE)
    vv <- 1
    if (cont) {
        nr.var <- 1
        if (cont && model.switch != 0) {
            vv <- ans.all$vv
            nr.var <- length(ans.all$vv)
        }
    }
    if (!cont) 
        nr.var <- 0
    if (model.switch == 0) {
        skip.b <- FALSE
        skip.0 <- FALSE
    }
    if (nr.covar > 10) 
        skip.0 <- T
    if (ans.all$DA.ans == 2) {
        skip.0 <- TRUE
        skip.a <- TRUE
        skip.ab <- TRUE
    }
    else skip.a <- FALSE
    dd.start <- 1
    if (model.switch != 0) {
        Vaa.start <- ans.all$Vaa.start
        aa.start <- ans.all$aa.start
        Vced.start <- ans.all$Vced.start
        ced.start <- ans.all$ced.start
        if (length(ans.all$xans) > 1) 
            RPF.start <- ans.all$RPF.start
        if (dtype == 3) {
            ced.start <- NA
            Vced.start <- rep(NA, length(Vced.start))
        }
        cc.start <- ans.all$cc.start
        if (cont) 
            var.start <- ans.all$var.start
        ths.start <- ans.all$ths.start
        if (vv == 2) 
            fct3 <- covariate
        else fct3 <- rep(1, length(x))
        if (cont) 
            nr.var <- max(fct3)
        if (dtype == 6) 
            alfa.start <- ans.all$alfa.start
        aic.m1 <- ans.all$aic.m1
        row.sel.m1 <- ans.all$row.sel.m1
        skip.0 <- ans.all$skip.0
        skip.b <- ans.all$skip.b
    }
    ii <- 0
    if (model.switch == 0) 
        title.tmp <- "nested exponential models"
    if (model.switch == 1) 
        title.tmp <- "nested Hill models"
    if (model.switch == 2) 
        title.tmp <- "nested inverse exponential models"
    if (model.switch == 3) 
        title.tmp <- "LN models"
    name.tmp <- ""
    if (!ans.all$WAPP && output && interrupt && cont)
        if (display_plots) {
            f.graph.window(12, WAPP = ans.all$WAPP, title = title.tmp, 
                name.wapp = name.tmp, plotprefix = ans.all$plotprefix)
        }
    tb <- "\t"
    model.txt <- character()
    Vconverged <- numeric()
    Vloglik <- numeric()
    Vaic <- numeric()
    Vnpar <- integer()
    if (cont) 
        MLE.all <- matrix(NA, 18, 3 * nr.covar + 4)
    if (dtype %in% 2:3) 
        MLE.all <- matrix(NA, 18, 2 * nr.covar + 2 + nth + 1)
    if (dtype == 4) 
        MLE.all <- matrix(NA, 18, 2 * nr.covar + 2 + 1 + 1)
    if (dtype == 6) 
        MLE.all <- matrix(NA, 18, 1 + 2 * nr.covar + 2 + 1 + 
            1)
    if (length(xans) > 1) {
        MLE.all <- cbind(MLE.all, rep(NA, 18))
        RPF.tmp <- matrix(NA, 18, ans.all$nr.dosecol - 1)
        MLE.all <- cbind(MLE.all, RPF.tmp)
    }
    ans.all.fit <- ans.all
    if (ans.all$displ.no == 0) 
        ans.all.fit$displ.fact <- covariate
    if (output) {
        if (model.switch == 0) 
            cat("\nANALYSIS WITH EXPONENTIAL MODELS\n")
        if (model.switch == 1) 
            cat(" \nANALYSIS WITH HILL MODELS\n")
        if (model.switch == 2) 
            cat(" \nANALYSIS WITH INVERSE EXPONENTIAL MODELS")
        if (model.switch == 3) 
            cat(" \nANALYSIS WITH LOGNORMAL DR MODELS")
        cat("\nmodel", tb, "converged", tb, "npar", tb, "loglik", 
            tb, "aic")
    }
    if (model.switch == 0) {
        fct3 <- rep(1, length(x))
        if (full.ans == 1) {
            ii <- ii + 1
            model.txt[ii] <- "full model"
            ans.all.fit$model.ans <- 11
            if (dtype %in% c(4, 6)) {
                ans.all.fit$model.type <- 1
                ans.all.fit$model.ans <- 14
            }
            ans.all.fit$modelname <- model.txt[ii]
            ans.all.fit$fct1 <- covariate
            ans.all.fit$fct2 <- covariate
            ans.all.fit$fct3 <- fct3
            if (dtype != 6) 
                ans.all.fit <- f.qfit(ans.all.fit, plot.type, 
                  output = output, display_plots = display_plots)
            else {
                ans.all.fit <- f.dtype6.mn(ans.all.fit)
                ans.all$full6.done <- TRUE
                ans.all$alfa.mle <- ans.all.fit$alfa.mle
                ans.all$Pvalue.alfa <- ans.all.fit$Pvalue.alfa
                ans.all$x.full <- ans.all.fit$x.full
                ans.all$pi.full <- ans.all.fit$pi.full
                ans.all$fct1.full <- ans.all.fit$fct1.full
                ans.all$fct2.full <- ans.all.fit$fct2.full
            }
            Vnpar[ii] <- ans.all.fit$npar.aic
            Vloglik[ii] <- round(ans.all.fit$loglik, 2)
            Vaic[ii] <- ans.all.fit$aic
            Vconverged[ii] <- ans.all.fit$converged
            MLE <- ans.all.fit$MLE
            MLE.all[ii, ] <- signif(MLE[1:length(MLE.all[1, ])], 
                7)
            if (output) 
                cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
                  tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            if (cont) 
                var.start <- MLE[1]
            if (dtype == 6) {
                alfa.start <- MLE[1]
                ans.all$alfa.start <- alfa.start
                alfa <- MLE[1]
            }
            results.full <- f.model.specific.results(ans.all.fit)
            if (ans.all$const.var == TRUE) 
                vv <- 1
            else if (cont && covar.no > 0) {
                ii <- ii + 1
                model.txt[ii] <- paste("full-", "v", sep = "")
                ans.all.fit$model.ans <- 11
                ans.all.fit$modelname <- model.txt[ii]
                fct3 <- covariate
                ans.all.fit$fct3 <- fct3
                if (cont) 
                  nr.var <- max(fct3)
                ans.all.fit$par.start <- c(rep(var.start, nr.var), 
                  ans.all.fit$par.start[-1])
                ans.all.fit <- f.qfit(ans.all.fit, plot.type = 0, 
                  output = output, display_plots = display_plots)
                Vnpar[ii] <- ans.all.fit$npar.aic
                Vloglik[ii] <- round(ans.all.fit$loglik, 2)
                Vaic[ii] <- ans.all.fit$aic
                Vconverged[ii] <- ans.all.fit$converged
                MLE <- ans.all$MLE
                if (dtype == 6) 
                  alfa <- MLE[1]
                if (Vaic[ii] > Vaic[ii - 1]) {
                  vv <- 1
                  fct3 <- rep(1, length(x))
                  nr.var <- 1
                }
                else {
                  vv <- 2
                  results.full <- f.model.specific.results(ans.all.fit)
                }
                MLE <- ans.all.fit$MLE
                MLE.all[ii, ] <- signif(MLE[1:length(MLE.all[1, 
                  ])], 7)
                if (output) 
                  cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
                    tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            }
            ans.all$FULL.model <- results.full
        }
        if (full.ans == 2) {
            ii <- 1
            model.txt[1] <- "full"
            if (cont) 
                nr.var <- 1
            if (dtype == 6) {
                Vnpar[ii] <- ans.all$npar
                Vloglik[ii] <- round(ans.all$loglik, 2)
                Vaic[ii] <- ans.all$aic
                Vconverged[ii] <- ans.all$converged
                if (output) 
                  cat(paste("\n", model.txt[ii], tb, tb, Vconverged[ii], 
                    tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
                ans.all$loglik.full <- Vloglik[ii]
                ans.all$aic.full <- Vaic[ii]
                ans.all$npar.aic.full <- Vnpar[ii]
            }
            else {
                Vloglik[1] <- NA
                Vaic[1] <- NA
                Vnpar[1] <- NA
                MLE.all[1, 1:4] <- 0
                if (output) 
                  cat(paste("\n", model.txt[ii], tb, "not evaluated"))
            }
        }
        if (!cont) 
            ans.all.fit$model.type <- 2
        if (length(xans) == 1) {
            ii <- ii + 1
            row.sel.m1 <- ii
            model.txt[ii] <- "null model"
            if (cont) 
                nr.var <- max(fct3)
            if (nr.var > 1) {
                model.txt[ii] <- paste(model.txt[ii], "v", sep = "-")
                ans.all.fit$nr.var <- nr.var
            }
            model.ans <- 1
            ans.all.fit$model.ans <- model.ans
            ans.all.fit$modelname <- model.txt[ii]
            fct1 <- rep(1, length(x))
            fct2 <- rep(1, length(x))
            ans.all.fit$fct1 <- fct1
            ans.all.fit$fct2 <- fct2
            ans.all.fit$fct3 <- fct3
            ans.all.fit$nr.aa <- 1
            ans.all.fit$nr.bb <- 1
            ans.all.fit$par.start <- NA
            if (full.ans == 2 || !cont) 
                var.start <- NA
            if (cont && full.ans == 1) 
                ans.all.fit$par.start <- c(rep(var.start, nr.var), 
                  exp(mean(log(y))))
            ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output, display_plots = display_plots)
            if (is.na(ans.all.fit$loglik)) {
                if (output) 
                  cat("\nATTENTION:  model 1 did not result in finite log-likelihood valu\n\n")
                ans.all.fit$loglik <- NA
            }
            Vnpar[ii] <- ans.all.fit$npar.aic
            Vloglik[ii] <- round(ans.all.fit$loglik, 2)
            Vaic[ii] <- ans.all.fit$aic
            Vconverged[ii] <- ans.all.fit$converged
            MLE <- ans.all.fit$MLE
            MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
            if (output) 
                cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
                  tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            if (cont && full.ans == 2) 
                var.start <- MLE[1]
            aa.start <- MLE[nr.var + 1]
            if (!cont) {
                npar.tmp <- length(MLE)
                ths.start <- MLE[(alfa.length + 1 + 1):(npar.tmp)]
            }
            aic.m1 <- ans.all.fit$aic
            results.m1 <- f.model.specific.results(ans.all.fit)
            if (!skip.a) 
                if (covar.no > 0) {
                  ii <- ii + 1
                  model.txt[ii] <- "null model-a"
                  if (max(fct3) > 1) 
                    model.txt[ii] <- paste(model.txt[ii], "v", 
                      sep = "-")
                  model.ans <- 1
                  ans.all.fit$model.ans <- model.ans
                  ans.all.fit$modelname <- model.txt[ii]
                  fct1 <- covariate
                  fct2 <- rep(1, length(x))
                  nr.aa <- max(fct1)
                  ans.all.fit$fct1 <- fct1
                  ans.all.fit$fct2 <- fct2
                  ans.all.fit$fct3 <- fct3
                  ans.all.fit$nr.aa <- nr.aa
                  ans.all.fit$nr.bb <- 1
                  if (cont && full.ans == 2) 
                    var.start <- NA
                  if (cont) 
                    par.start <- c(rep(var.start, nr.var), rep(aa.start, 
                      nr.aa))
                  if (!cont) 
                    par.start <- c(rep(aa.start, nr.aa), ths.start)
                  if (dtype == 6) 
                    par.start <- c(alfa.start, par.start)
                  ans.all.fit$par.start <- par.start
                  ans.all.fit <- f.qfit(ans.all.fit, plot.type, 
                    output = output, display_plots = display_plots)
                  Vnpar[ii] <- ans.all.fit$npar.aic
                  Vloglik[ii] <- round(ans.all.fit$loglik, 2)
                  Vaic[ii] <- ans.all.fit$aic
                  Vconverged[ii] <- ans.all.fit$converged
                  MLE <- ans.all.fit$MLE
                  MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
                  if (is.na(Vloglik[ii])) {
                    if (output) 
                      cat("\nATTENTION:  model 1 did not result in finite log-likelihood value \n\n")
                    Vloglik[ii] <- NA
                  }
                  if (output) 
                    cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
                      tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
                  if (Vaic[ii] < Vaic[ii - 1]) {
                    results.m1 <- f.model.specific.results(ans.all.fit)
                    row.sel.m1 <- ii
                    aic.m1 <- ans.all.fit$aic
                  }
                  Vaa.start <- MLE[(nr.var + 1):(nr.var + nr.aa)]
                }
            ans.all$NULL.model <- results.m1
        }
        if (length(ans.all$xans) > 1 || skip.a) {
            aa.start <- NA
            Vaa.start <- NA
            if (cont) 
                var.start <- NA
            if (output) 
                cat(paste("\n null", tb, "not evaluated"))
            aic.m1 <- NA
            row.sel.m1 <- NA
            RPF.start <- NA
        }
        ans.all.fit$model.ans <- 2
        ans.all.fit$modelname <- "m2"
        if (max(covariate) > 1) 
            ans.all.fit$modelname <- "m2-ab"
        if (max(fct3) > 1) 
            ans.all.fit$modelname <- paste(ans.all.fit$modelname, 
                "v", sep = "")
        fct1 <- covariate
        fct2 <- covariate
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        if (cont) 
            nr.var <- max(fct3)
        ans.all.fit$fct1 <- fct1
        ans.all.fit$fct2 <- fct2
        if (vv == 1) 
            ans.all.fit$fct3 <- rep(1, length(x))
        if (vv == 2) 
            ans.all.fit$fct3 <- covariate
        ans.all.fit$nr.aa <- nr.aa
        ans.all.fit$nr.bb <- nr.bb
        if (cont) 
            ans.all.fit$nr.var <- max(ans.all.fit$fct3)
        if (cont) 
            par.start.tmp <- f.start.con(ans.all.fit, tmp.quick = T)$par.start
        if (!cont) 
            par.start.tmp <- f.start.cat(ans.all.fit, tmp.quick = T)$par.start
        if (covar.no == 0) 
            aa.start.tmp <- aa.start
        if (covar.no > 0) 
            aa.start.tmp <- Vaa.start
        bb.start <- par.start.tmp[(alfa.length + nr.var + nr.aa + 
            1):(alfa.length + nr.var + nr.aa + nr.bb)]
        if (cont) 
            ans.all.fit$par.start <- c(rep(var.start, nr.var), 
                aa.start.tmp, bb.start)
        else ans.all.fit$par.start <- NA
        ans.all.fit <- f.qfit(ans.all.fit, plot.type = 0, output = output, display_plots = display_plots)
        MLE <- ans.all.fit$MLE
        if (0) {
            Vnpar[ii] <- ans.all.fit$npar.aic
            Vloglik[ii] <- round(ans.all.fit$loglik, 2)
            Vaic[ii] <- ans.all.fit$aic
            Vconverged[ii] <- ans.all.fit$converged
            if (vv == 2) 
                cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                  tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            else cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
        }
        if (output) 
            if (length(xans) > 1) 
                cat(paste(" \n", "model 2", tb, ans.all.fit$loglik))
        Vaa.start <- MLE[(alfa.length + nr.var + 1):(alfa.length + 
            nr.var + nr.aa)]
        aa.start <- exp(mean(log(Vaa.start)))
        Vbb.start <- MLE[(alfa.length + nr.var + max(fct1) + 
            1):(alfa.length + nr.var + max(fct1) + max(fct2))]
        if (!cont) {
            npar.tmp <- length(MLE)
            ths.start <- MLE[(alfa.length + nr.aa + nr.bb + 1):npar.tmp]
        }
        if (cont && full.ans == 2) 
            var.start <- MLE[alfa.length + 1]
        if (length(xans) > 1) {
            RPF.start <- MLE[(length(MLE) - ans.all$nr.dosecol + 
                2):length(MLE)]
            ans.all$RPF.start <- RPF.start
        }
        if (output) 
            if (dtype != 6 & alfa.length != 0) {
                print("f.select.m5.con")
                cat(" \n\nalfa.length not equal to zero \n\n")
                f.press.key.to.continue()
            }
        if (cont) {
            if (ans.all$MA.running) 
                increase <- ans.all$increase
            else {
                if (abs(sum(sign(Vbb.start))) != length(Vbb.start)) 
                  cat(" \nSlopes have opposite signs")
                if (mean(Vbb.start) > 0) {
                  increase <- 1
                  CES <- abs(CES)
                }
                else {
                  increase <- -1
                  CES <- -abs(CES)
                }
                if (CES <= -1) {
                  CES.old <- CES
                  FC <- abs(CES) + 1
                  CES <- 1/FC - 1
                  cat(" \nATTENTION: value of CES (", CES.old, 
                    ") has been changed to", CES, "\n               as original value is not possible for decreasing dose-responses \n\n")
                }
                ans.all.fit$increase <- increase
                ans.all.fit$CES <- CES
                ans.all$increase <- increase
                ans.all$CES <- CES
            }
        }
        if (cont) {
            regr.par.matr <- ans.all.fit$regr.par.matr
            bb.tmp <- regr.par.matr[, 2]
            if (ans.all$increase == 1) {
                bb.tmp[bb.tmp < 0] <- min(bb.tmp[bb.tmp > 0])
                regr.par.matr[, 2] <- bb.tmp
            }
            if (ans.all$increase == -1) {
                bb.tmp[bb.tmp > 0] <- max(bb.tmp[bb.tmp < 0])
                regr.par.matr[, 2] <- bb.tmp
            }
            ans.all.ced <- ans.all
            ans.all.ced$model.ans <- 2
            ans.all.ced$regr.par.matr <- regr.par.matr
            ans.all.ced$nr.aa <- nr.aa
            ans.all.ced$nr.bb <- nr.bb
            ans.all.ced$nr.var <- 1
            ced.lst <- f.ced.con(ans.all.ced)
            Vced.start <- ced.lst$CED
            ced.start <- exp(mean(log(Vced.start)))
        }
        if (!cont) {
            ced.lst <- f.ced.cat(ans.all.fit)
            Vced.start <- ced.lst$CED.matr[, 1]
            if (all(is.finite(Vced.start))) {
                ced.start <- exp(mean(log(Vced.start)))
                Vced.start[Vced.start == Inf] <- Vced.start
            }
            else {
                Vced.start <- rep(NA, length(Vced.start))
                ced.start <- NA
            }
        }
        if (ans.all$increase == -1) 
            cc.start <- 0.1
        if (ans.all$increase == 1) 
            cc.start <- 10
        if (!cont) 
            cc.start <- 0.1
        ces.cat <- ans.all$CES.cat
        ces.cat.tmp <- 1
        if (dtype == 3) {
            for (ii in 1:ans.all.fit$nr.gr) {
                pi.matr <- f.expect.cat(model.type = 2, model.ans = 2, 
                  x = 0, ans.all.fit$regr.par.matr[ii, ], ans.all.fit$th.par, 
                  ans.all.fit$sig.par, CES = 0, ces.ans = 1, 
                  dtype = 3, twice = TRUE, cc.inf = cc.inf)
                dum <- pi.matr < 0.5
                ces.cat.tmp <- min(order(dum)[dum == TRUE])
                ces.cat <- max(ces.cat.tmp, ces.cat)
            }
            if (ces.cat > ans.all$CES.cat) {
                cat("\nATTENTION: Severity category for CED was set to", 
                  ces.cat, "due to properties of the data \n\n")
                ans.all.fit$CES.cat <- ces.cat
                ans.all$CES.cat <- ces.cat
            }
        }
    }
    ans.all.fit$increase <- ans.all$increase
    if (model.switch == 0) 
        nr.first.models <- ii
    else nr.first.models <- 0
    ii.m3 <- 0
    m3x <- numeric()
    Vaic.m3 <- numeric()
    model.ans <- 13 + (model.switch == 1) * 10 + (model.switch == 
        2) * 38 + (model.switch == 3) * 40
    if (!skip.0) {
        ii <- ii + 1
        ii.m3 <- ii.m3 + 1
        if (model.switch == 0) 
            model.txt[ii] <- "Expon. m3-"
        if (model.switch == 1) 
            model.txt[ii] <- "Hill m3-"
        if (model.switch == 2) 
            model.txt[ii] <- "Inv.Expon. m3-"
        if (model.switch == 3) 
            model.txt[ii] <- "LN m3-"
        fct1 <- rep(1, length(x))
        fct2 <- rep(1, length(x))
        if (vv == 2) 
            fct3 <- covariate
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        if (cont) 
            nr.var <- max(fct3)
        if (max(fct3) > 1) 
            model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
        ans.all.fit$model.ans <- model.ans
        ans.all.fit$modelname <- model.txt[ii]
        ans.all.fit$fct1 <- fct1
        ans.all.fit$fct2 <- fct2
        ans.all.fit$fct3 <- fct3
        ans.all.fit$nr.aa <- nr.aa
        ans.all.fit$nr.bb <- nr.bb
        ans.all.fit$nr.var <- nr.var
        if (cont) 
            par.start <- c(rep(var.start, nr.var), aa.start, 
                ced.start, 1)
        if (!cont) 
            par.start <- c(aa.start, ced.start, 1, ths.start)
        if (dtype == 6) 
            par.start <- c(alfa.start, par.start)
        if (length(xans) > 1) 
            par.start <- c(par.start, RPF.start)
        ans.all.fit$par.start <- par.start
        ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output, display_plots = display_plots)
        Vnpar[ii] <- ans.all.fit$npar.aic
        Vloglik[ii] <- round(ans.all.fit$loglik, 2)
        Vaic[ii] <- ans.all.fit$aic
        Vconverged[ii] <- ans.all.fit$converged
        MLE <- ans.all.fit$MLE
        MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
        Vaic.m3[ii.m3] <- Vaic[ii]
        m3x[ii.m3] <- "0"
        if (output) {
            cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
        }
        aa.start <- MLE[alfa.length + nr.var + 1]
        if (cont && full.ans == 2) 
            var.start <- MLE[1]
        if (length(xans) > 1) 
            RPF.start <- MLE[(length(MLE) - ans.all$nr.dosecol + 
                2):length(MLE)]
        results.m3.0 <- f.model.specific.results(ans.all.fit)
    }
    if (covar.no > 0) {
        if (!skip.a) {
            ii <- ii + 1
            ii.m3 <- ii.m3 + 1
            if (model.switch == 0) 
                model.txt[ii] <- "Expon. m3-a"
            if (model.switch == 1) 
                model.txt[ii] <- "Hill m3-a"
            if (model.switch == 2) 
                model.txt[ii] <- "Inv.Expon. m3-a"
            if (model.switch == 3) 
                model.txt[ii] <- "LN m3-a"
            if (max(fct3) > 1) 
                model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
            ans.all.fit$model.ans <- model.ans
            ans.all.fit$modelname <- model.txt[ii]
            fct1 <- covariate
            ans.all.fit$fct1 <- fct1
            ans.all.fit$fct2 <- rep(1, length(x))
            ans.all.fit$fct3 <- fct3
            nr.aa <- max(fct1)
            nr.bb <- 1
            ans.all.fit$nr.aa <- nr.aa
            ans.all.fit$nr.bb <- nr.bb
            ans.all.fit$nr.var <- nr.var
            if (cont) 
                par.start <- c(rep(var.start, nr.var), Vaa.start, 
                  ced.start, 1)
            if (!cont) 
                par.start <- c(aa.start, Vced.start, 1, ths.start)
            if (dtype == 6) 
                par.start <- c(alfa.start, par.start)
            if (length(xans) > 1) 
                par.start <- c(par.start, RPF.start)
            ans.all.fit$par.start <- par.start
            ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output, display_plots = display_plots)
            Vnpar[ii] <- ans.all.fit$npar.aic
            Vloglik[ii] <- round(ans.all.fit$loglik, 2)
            Vaic[ii] <- ans.all.fit$aic
            Vconverged[ii] <- ans.all.fit$converged
            MLE <- ans.all.fit$MLE
            MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
            Vaic.m3[ii.m3] <- Vaic[ii]
            m3x[ii.m3] <- "a"
            if (output) 
                cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                  tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            Vaa.start <- MLE[(alfa.length + nr.var + 1):(alfa.length + 
                nr.var + max(fct1))]
            results.m3.a <- f.model.specific.results(ans.all.fit)
            if (model.switch == 0) 
                if (Vaic[ii - 1] - Vaic[ii] > 10) {
                  skip.b <- T
                  skip.0 <- T
                }
                else {
                  skip.b <- F
                }
            if (length(xans) > 1) {
                skip.b <- TRUE
                skip.ab <- TRUE
            }
            else skip.ab <- FALSE
        }
        if (!skip.b) {
            ii <- ii + 1
            ii.m3 <- ii.m3 + 1
            if (model.switch == 0) 
                model.txt[ii] <- "Expon. m3-b"
            if (model.switch == 1) 
                model.txt[ii] <- "Hill m3-b"
            if (model.switch == 2) 
                model.txt[ii] <- "Inv.Expon. m3-b"
            if (model.switch == 3) 
                model.txt[ii] <- "LN m3-b"
            if (max(fct3) > 1) 
                model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
            ans.all.fit$model.ans <- model.ans
            ans.all.fit$modelname <- model.txt[ii]
            fct1 <- rep(1, length(x))
            fct2 <- covariate
            ans.all.fit$fct1 <- fct1
            ans.all.fit$fct2 <- fct2
            ans.all.fit$fct3 <- fct3
            nr.aa <- 1
            nr.bb <- max(fct2)
            ans.all.fit$nr.aa <- nr.aa
            ans.all.fit$nr.bb <- nr.bb
            ans.all.fit$nr.var <- nr.var
            if (cont) 
                par.start <- c(rep(var.start, nr.var), aa.start, 
                  Vced.start, 1)
            if (!cont && dtype != 3) 
                par.start <- c(aa.start, Vced.start, 1, ths.start)
            if (dtype == 3) 
                par.start <- NA
            if (dtype == 6) 
                par.start <- c(alfa.start, par.start)
            ans.all.fit$par.start <- par.start
            ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output, display_plots = display_plots)
            Vnpar[ii] <- ans.all.fit$npar.aic
            Vloglik[ii] <- round(ans.all.fit$loglik, 2)
            Vaic[ii] <- ans.all.fit$aic
            Vconverged[ii] <- ans.all.fit$converged
            MLE <- ans.all.fit$MLE
            MLE.2.b <- MLE
            MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
            Vaic.m3[ii.m3] <- Vaic[ii]
            m3x[ii.m3] <- "b"
            if (output) 
                cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                  tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            aa.start.mb <- MLE[alfa.length + nr.var + 1]
            results.m3.b <- f.model.specific.results(ans.all.fit)
        }
        if (!skip.ab) {
            ii <- ii + 1
            ii.m3 <- ii.m3 + 1
            if (model.switch == 0) 
                model.txt[ii] <- "Expon. m3-ab"
            if (model.switch == 1) 
                model.txt[ii] <- "Hill m3-ab"
            if (model.switch == 2) 
                model.txt[ii] <- "Inv.Expon. m3-ab"
            if (model.switch == 3) 
                model.txt[ii] <- "LN m3-ab"
            if (max(fct3) > 1) 
                model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
            ans.all.fit$modelname <- model.txt[ii]
            fct1 <- covariate
            fct2 <- covariate
            ans.all.fit$fct1 <- fct1
            ans.all.fit$fct2 <- fct2
            ans.all.fit$fct3 <- fct3
            nr.aa <- max(fct1)
            nr.bb <- max(fct2)
            if (cont) 
                nr.var <- max(fct3)
            ans.all.fit$nr.aa <- nr.aa
            ans.all.fit$nr.bb <- nr.bb
            ans.all.fit$nr.var <- nr.var
            if (length(Vced.start) == 1) 
                Vced.start <- rep(Vced.start, nr.bb)
            if (cont) 
                par.start <- c(rep(var.start, nr.var), rep(aa.start, 
                  nr.aa), Vced.start, 1)
            if (!cont && dtype != 3) 
                par.start <- c(rep(aa.start, nr.aa), Vced.start, 
                  1, ths.start)
            if (dtype == 3) 
                par.start <- NA
            if (dtype == 6) 
                par.start <- c(alfa.start, par.start)
            ans.all.fit$par.start <- par.start
            ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output, display_plots = display_plots)
            Vnpar[ii] <- ans.all.fit$npar.aic
            Vloglik[ii] <- round(ans.all.fit$loglik, 2)
            Vaic[ii] <- ans.all.fit$aic
            Vconverged[ii] <- ans.all.fit$converged
            MLE <- ans.all.fit$MLE
            MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
            Vaic.m3[ii.m3] <- Vaic[ii]
            m3x[ii.m3] <- "ab"
            if (output) {
                if (vv == 2) 
                  cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                    tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
                else cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                  tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            }
            Vaa.start <- MLE[(alfa.length + nr.var + 1):(alfa.length + 
                nr.var + max(fct1))]
            Vced.start <- MLE[(alfa.length + nr.var + max(fct1) + 
                1):(alfa.length + nr.var + max(fct1) + max(fct2))]
            if (!cont) {
                npar.tmp <- length(MLE)
                ths.start <- MLE[(alfa.length + nr.aa + nr.bb + 
                  2):(npar.tmp)]
            }
            if (dtype != 6 & alfa.length != 0) {
                print("f.select.m5.con")
                cat(" \nalfa.length not equal to zero \n\n")
                f.press.key.to.continue()
            }
            results.m3.ab <- f.model.specific.results(ans.all.fit)
        }
    }
    rank.m3 <- order(Vaic.m3)[1]
    row.sel.m3 <- nr.first.models + rank.m3
    if (quick.ans == 4) 
        if (covar.no > 0) {
            if (output) {
                cat(" \n--------------------------------------------------- \n")
                cat("Best model 3 with covariates is:", model.txt[row.sel.m3], 
                  " \n\n")
            }
            if (ii.m3 > 1) {
                rank.m3.2 <- order(Vaic.m3)[2]
                if (output && Vaic.m3[rank.m3.2] - Vaic.m3[rank.m3] <= 
                  2) 
                  cat("However", model.txt[nr.first.models + 
                    rank.m3.2], "is a reasonable model as well \n")
            }
            if (ii.m3 > 2) {
                rank.m3.3 <- order(Vaic.m3)[3]
                if (output && Vaic.m3[rank.m3.3] - Vaic.m3[rank.m3] <= 
                  2) 
                  cat("and so is", model.txt[nr.first.models + 
                    rank.m3.3], " \n")
            }
            if (ii.m3 > 3) {
                rank.m3.4 <- order(Vaic.m3)[3]
                if (output && Vaic.m3[rank.m3.4] - Vaic.m3[rank.m3] <= 
                  2) 
                  cat("and so is", model.txt[nr.first.models + 
                    rank.m3.4], " \n\n")
            }
        }
    if (quick.ans != 4) {
        ii.m5 <- 0
        m5x <- numeric()
        Vaic.m5 <- numeric()
        nr.first.models.m3or5 <- nr.first.models
        nr.first.models <- ii
        model.ans <- 15 + (model.switch == 1) * 10 + (model.switch == 
            2) * 37 + (model.switch == 3) * 39
        if (!skip.0) {
            ii <- ii + 1
            ii.m5 <- ii.m5 + 1
            if (model.switch == 0) 
                model.txt[ii] <- "Expon. m5-"
            if (model.switch == 1) 
                model.txt[ii] <- "Hill m5-"
            if (model.switch == 2) 
                model.txt[ii] <- "Inv.Expon. m5-"
            if (model.switch == 3) 
                model.txt[ii] <- "LN m5-"
            fct1 <- rep(1, length(x))
            fct2 <- rep(1, length(x))
            if (max(fct3) > 1) 
                model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
            ans.all.fit$model.ans <- model.ans
            ans.all.fit$modelname <- model.txt[ii]
            ans.all.fit$fct1 <- fct1
            ans.all.fit$fct2 <- fct2
            ans.all.fit$fct3 <- fct3
            ans.all.fit$nr.aa <- max(fct1)
            ans.all.fit$nr.bb <- max(fct2)
            ans.all.fit$nr.var <- max(fct3)
            if (cont) 
                nr.var <- max(fct3)
            ans.all.fit$nr.var <- nr.var
            if (0) {
                if (cont) 
                  par.start <- c(rep(var.start, nr.var), aa.start, 
                    ced.start, cc.start, dd.start)
                if (!cont && dtype != 3) 
                  par.start <- c(aa.start, ced.start, cc.start, 
                    dd.start, ths.start)
                if (dtype == 3) 
                  par.start <- NA
                if (dtype == 6) 
                  par.start <- c(alfa.start, par.start)
                if (length(xans) > 1) 
                  par.start <- c(par.start, RPF.start)
                ans.all.fit$par.start <- par.start
            }
            ans.all.fit$par.start <- NA
            ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output, display_plots = display_plots)
            Vnpar[ii] <- ans.all.fit$npar.aic
            Vloglik[ii] <- round(ans.all.fit$loglik, 2)
            Vaic[ii] <- ans.all.fit$aic
            Vconverged[ii] <- ans.all.fit$converged
            MLE <- ans.all.fit$MLE
            MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
            Vaic.m5[ii.m5] <- Vaic[ii]
            m5x[ii.m5] <- "0"
            if (output) 
                cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                  tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            results.m5.0 <- f.model.specific.results(ans.all.fit)
        }
        if (covar.no > 0) {
            if (!skip.a) {
                ii <- ii + 1
                ii.m5 <- ii.m5 + 1
                model.ans <- 15 + (model.switch == 1) * 10 + 
                  (model.switch == 2) * 37 + (model.switch == 
                  3) * 39
                model.txt[ii] <- "Expon. m5-a"
                if (model.switch == 1) 
                  model.txt[ii] <- "Hill m5-a"
                if (model.switch == 2) 
                  model.txt[ii] <- "Inv.Expon. m5-a"
                if (model.switch == 3) 
                  model.txt[ii] <- "LN m5-a"
                if (max(fct3) > 1) 
                  model.txt[ii] <- paste(model.txt[ii], "v", 
                    sep = "")
                ans.all.fit$model.ans <- model.ans
                ans.all.fit$modelname <- model.txt[ii]
                fct1 <- covariate
                fct2 <- rep(1, length(x))
                ans.all.fit$fct1 <- fct1
                ans.all.fit$fct2 <- fct2
                ans.all.fit$fct3 <- fct3
                nr.aa <- max(fct1)
                nr.bb <- max(fct2)
                ans.all.fit$nr.aa <- nr.aa
                ans.all.fit$nr.bb <- nr.bb
                ans.all.fit$nr.var <- nr.var
                if (cont) 
                  par.start <- c(rep(var.start, nr.var), Vaa.start, 
                    ced.start, cc.start, dd.start)
                if (!cont && dtype != 3) 
                  par.start <- c(Vaa.start, ced.start, cc.start, 
                    dd.start, ths.start)
                if (dtype == 3) 
                  par.start <- NA
                if (length(xans) > 1) 
                  par.start <- c(par.start, RPF.start)
                if (dtype == 6) 
                  par.start <- c(alfa.start, par.start)
                ans.all.fit$par.start <- par.start
                ans.all.fit <- f.qfit(ans.all.fit, plot.type, 
                  output = output, display_plots = display_plots)
                Vnpar[ii] <- ans.all.fit$npar.aic
                Vloglik[ii] <- round(ans.all.fit$loglik, 2)
                Vaic[ii] <- ans.all.fit$aic
                Vconverged[ii] <- ans.all.fit$converged
                MLE <- ans.all.fit$MLE
                MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
                Vaic.m5[ii.m5] <- Vaic[ii]
                m5x[ii.m5] <- "a"
                if (output) 
                  cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                    tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
                Vaa.start <- MLE[(alfa.length + nr.var + 1):(alfa.length + 
                  nr.var + max(fct1))]
                results.m5.a <- f.model.specific.results(ans.all.fit)
            }
            if (!skip.b) {
                ii <- ii + 1
                ii.m5 <- ii.m5 + 1
                if (model.switch == 0) 
                  model.txt[ii] <- "Expon. m5-b"
                if (model.switch == 1) 
                  model.txt[ii] <- "Hill m5-b"
                if (model.switch == 2) 
                  model.txt[ii] <- "Inv.Expon. m5-b"
                if (model.switch == 3) 
                  model.txt[ii] <- "LN m5-b"
                if (max(fct3) > 1) 
                  model.txt[ii] <- paste(model.txt[ii], "v", 
                    sep = "")
                ans.all.fit$model.ans <- model.ans
                ans.all.fit$modelname <- model.txt[ii]
                fct1 <- rep(1, length(x))
                fct2 <- covariate
                ans.all.fit$fct1 <- fct1
                ans.all.fit$fct2 <- fct2
                ans.all.fit$fct3 <- fct3
                nr.aa <- 1
                nr.bb <- max(fct2)
                ans.all.fit$nr.aa <- nr.aa
                ans.all.fit$nr.bb <- nr.bb
                ans.all.fit$nr.var <- nr.var
                if (cont) 
                  par.start <- c(rep(var.start, nr.var), aa.start.mb, 
                    Vced.start, cc.start, dd.start)
                if (!cont && dtype != 3) 
                  par.start <- c(aa.start.mb, Vced.start, cc.start, 
                    dd.start, ths.start)
                if (dtype == 3) 
                  par.start <- NA
                if (dtype == 6) 
                  par.start <- c(alfa.start, par.start)
                ans.all.fit$par.start <- par.start
                ans.all.fit <- f.qfit(ans.all.fit, plot.type, 
                  output = output, display_plots = display_plots)
                Vnpar[ii] <- ans.all.fit$npar.aic
                Vloglik[ii] <- round(ans.all.fit$loglik, 2)
                Vaic[ii] <- ans.all.fit$aic
                Vconverged[ii] <- ans.all.fit$converged
                MLE <- ans.all.fit$MLE
                MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
                Vaic.m5[ii.m5] <- Vaic[ii]
                m5x[ii.m5] <- "b"
                if (output) 
                  cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                    tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
                results.m5.b <- f.model.specific.results(ans.all.fit)
            }
            if (!skip.ab) {
                ii <- ii + 1
                ii.m5 <- ii.m5 + 1
                if (model.switch == 0) 
                  model.txt[ii] <- "Expon. m5-ab"
                if (model.switch == 1) 
                  model.txt[ii] <- "Hill m5-ab"
                if (model.switch == 2) 
                  model.txt[ii] <- "Inv.Expon. m5-ab"
                if (model.switch == 3) 
                  model.txt[ii] <- "LN m5-ab"
                if (max(fct3) > 1) 
                  model.txt[ii] <- paste(model.txt[ii], "v", 
                    sep = "")
                ans.all.fit$model.ans <- model.ans
                ans.all.fit$modelname <- model.txt[ii]
                fct1 <- covariate
                fct2 <- covariate
                ans.all.fit$fct1 <- fct1
                ans.all.fit$fct2 <- fct2
                ans.all.fit$fct3 <- fct3
                nr.aa <- max(fct1)
                nr.bb <- max(fct2)
                ans.all.fit$nr.aa <- nr.aa
                ans.all.fit$nr.bb <- nr.bb
                ans.all.fit$nr.var <- nr.var
                if (cont) 
                  nr.var <- max(fct3)
                ans.all.fit$nr.var <- nr.var
                if (length(Vced.start) == 1) 
                  Vced.start <- rep(Vced.start, nr.bb)
                if (cont) 
                  par.start <- c(rep(var.start, nr.var), Vaa.start, 
                    Vced.start, cc.start, dd.start)
                if (!cont && dtype != 3) 
                  par.start <- c(Vaa.start, Vced.start, cc.start, 
                    dd.start, ths.start)
                if (dtype == 3) 
                  par.start <- NA
                if (dtype == 6) 
                  par.start <- c(alfa.start, par.start)
                ans.all.fit$par.start <- par.start
                ans.all.fit <- f.qfit(ans.all.fit, plot.type, 
                  output = output, display_plots = display_plots)
                Vnpar[ii] <- ans.all.fit$npar.aic
                Vloglik[ii] <- round(ans.all.fit$loglik, 2)
                Vaic[ii] <- ans.all.fit$aic
                Vconverged[ii] <- ans.all.fit$converged
                MLE <- ans.all.fit$MLE
                MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
                Vaic.m5[ii.m5] <- Vaic[ii]
                m5x[ii.m5] <- "ab"
                if (output) {
                  if (vv == 2) 
                    cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                      tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
                  else cat(paste(" \n", model.txt[ii], tb, Vconverged[ii], 
                    tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
                }
                results.m5.ab <- f.model.specific.results(ans.all.fit)
            }
        }
        rank.m5 <- order(Vaic.m5)[1]
        row.sel.m5 <- nr.first.models + rank.m5
        if (quick.ans == 5) 
            if (covar.no > 0) {
                if (output) {
                  cat(" \n--------------------------------------------------- \n\n")
                  cat("Best model 5 with covariates is:", model.txt[row.sel.m5], 
                    " \n\n")
                }
                if (ii.m5 > 1) {
                  rank.m5.2 <- order(Vaic.m5)[2]
                  if (output && Vaic.m5[rank.m5.2] - Vaic.m5[rank.m5] <= 
                    2) 
                    cat("However", model.txt[nr.first.models + 
                      rank.m5.2], "is a reasonable model as well \n\n")
                }
                if (ii.m5 > 2) {
                  rank.m5.3 <- order(Vaic.m5)[3]
                  if (output && Vaic.m5[rank.m5.3] - Vaic.m5[rank.m5] <= 
                    2) 
                    cat("and so is", model.txt[nr.first.models + 
                      rank.m5.3], " \n\n")
                }
                if (ii.m5 > 3) {
                  rank.m5.4 <- order(Vaic.m5)[3]
                  if (output && Vaic.m5[rank.m5.4] - Vaic.m5[rank.m5] <= 
                    2) 
                    cat("and so is", model.txt[nr.first.models + 
                      rank.m5.4], " \n\n")
                }
                if (output) 
                  cat("--------------------------------------------------- \n")
            }
    }
    if (quick.ans == 4) 
        Vaic.models <- Vaic.m3
    if (quick.ans == 5) 
        Vaic.models <- Vaic.m5
    if (quick.ans == 3) 
        Vaic.models <- c(Vaic.m3, Vaic.m5)
    model.sel <- 1000
    if (length(xans) == 1 && !skip.a) {
        if (all(Vaic.models - aic.m1 > -2)) {
            cat(" \nNo significant trend according to AIC criterion for model", 
                model.ans, " \n\n\n")
            model.sel <- 1
            row.sel <- row.sel.m1
        }
    }
    if (quick.ans == 3 && covar.no > 0 && model.sel > 1) {
        rank.first <- order(Vaic.models)[1]
        row.sel <- rank.first + nr.first.models.m3or5
        nr.first.models <- nr.first.models.m3or5
        if (covar.no > 0) {
            if (output) {
                cat("\n ---------------------------------------------------\n\n ")
                cat("Best model with covariates is:", model.txt[nr.first.models + 
                  rank.first], "\n\n ")
            }
            rank.second <- order(Vaic.models)[2]
            if (output && Vaic.models[rank.second] <= Vaic.models[rank.first] + 
                2) 
                cat("However", model.txt[nr.first.models + rank.second], 
                  "is a reasonable model as well\n\n ")
            if (length(Vaic.models) > 2) {
                rank.third <- order(Vaic.models)[3]
                if (output && Vaic.models[rank.third] <= Vaic.models[rank.first] + 
                  2) 
                  cat("and so is", model.txt[nr.first.models + 
                    rank.third], "\n\n ")
            }
            if (length(Vaic.models) > 3) {
                rank.fourth <- order(Vaic.models)[4]
                if (output && Vaic.models[rank.fourth] <= Vaic.models[rank.first] + 
                  2) 
                  cat("and so is", model.txt[nr.first.models + 
                    rank.fourth], "\n\n ")
            }
            if (length(Vaic.models) > 4) {
                rank.fifth <- order(Vaic.models)[5]
                if (output && Vaic.models[rank.fifth] <= Vaic.models[rank.first] + 
                  2) 
                  cat("and so is", model.txt[nr.first.models + 
                    rank.fifth], "\n\n ")
            }
        }
    }
    if (quick.ans == 4) {
        m3x.sel <- m3x[rank.m3]
        row.sel <- row.sel.m3
        if (m3x.sel == "0") 
            ans.all.spec <- results.m3.0
        if (m3x.sel == "a") 
            ans.all.spec <- results.m3.a
        if (m3x.sel == "b") 
            ans.all.spec <- results.m3.b
        if (m3x.sel == "ab") 
            ans.all.spec <- results.m3.ab
    }
    if (quick.ans == 5) {
        m5x.sel <- m5x[rank.m5]
        row.sel <- row.sel.m5
        if (m5x.sel == "0") 
            ans.all.spec <- results.m5.0
        if (m5x.sel == "a") 
            ans.all.spec <- results.m5.a
        if (m5x.sel == "b") 
            ans.all.spec <- results.m5.b
        if (m5x.sel == "ab") 
            ans.all.spec <- results.m5.ab
    }
    aic.tmp <- 0
    if (quick.ans == 3) {
        if (model.ans %in% c(4, 6)) 
            aic.tmp <- 5
        m5x.sel <- m5x[rank.m5]
        m3x.sel <- m3x[rank.m3]
        aic.m5.sel <- Vaic.m5[rank.m5]
        aic.m3.sel <- Vaic.m3[rank.m3]
        if (aic.m3.sel <= aic.m5.sel + aic.tmp) {
            row.sel <- row.sel.m3
            if (m3x.sel == "0") 
                ans.all.spec <- results.m3.0
            if (m3x.sel == "a") 
                ans.all.spec <- results.m3.a
            if (m3x.sel == "b") 
                ans.all.spec <- results.m3.b
            if (m3x.sel == "ab") 
                ans.all.spec <- results.m3.ab
        }
        if (aic.m5.sel + aic.tmp < aic.m3.sel) {
            row.sel <- row.sel.m5
            if (m5x.sel == "0") 
                ans.all.spec <- results.m5.0
            if (m5x.sel == "a") 
                ans.all.spec <- results.m5.a
            if (m5x.sel == "b") 
                ans.all.spec <- results.m5.b
            if (m5x.sel == "ab") 
                ans.all.spec <- results.m5.ab
        }
    }
    if (model.sel == 1) 
        ans.all.spec$trend <- FALSE
    else ans.all.spec$trend <- TRUE
    ans.all.tmp <- f.move.sublist(ans.all, ans.all.spec)
    if (!ans.all$const.var) 
        if (cont && full.ans == 2 && max(fct3) == 1 && covar.no > 
            0) {
            ii <- ii + 1
            model.txt[ii] <- paste(model.txt[row.sel], "v", sep = "")
            ans.all.tmp$fct3 <- covariate
            nr.var <- max(covariate)
            ans.all.tmp$nr.var <- nr.var
            par.start <- c(rep(var.start, nr.var), ans.all.tmp$MLE[-1])
            ans.all.tmp$par.start <- par.start
            if (length(Vloglik) > 9) 
                plot.type.tmp <- 0
            else plot.type.tmp <- plot.type
            ans.all.tmp <- f.qfit(ans.all.tmp, plot.type.tmp, 
                output = output, display_plots = display_plots)
            Vnpar[ii] <- ans.all.tmp$npar
            Vloglik[ii] <- round(ans.all.tmp$loglik, 2)
            Vaic[ii] <- ans.all.tmp$aic
            Vconverged[ii] <- ans.all.tmp$converged
            MLE <- ans.all.tmp$MLE
            MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
            if (output) 
                cat(paste("\n ", model.txt[ii], tb, Vconverged[ii], 
                  tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
            if (Vaic[ii] >= Vaic[row.sel]) {
                fct3 <- rep(1, length(x))
                nr.var <- 1
            }
            else {
                ans.all.spec <- f.model.specific.results(ans.all.tmp)
                row.sel <- ii
                vv <- 2
                if (model.sel == 1) {
                  cat("ATTENTION: Without covariate-dependent variances no trend was assessed (based on AIC).\n                   You might consider to repeat the analysis including the full model to check if this changes\n\n")
                }
            }
        }
    if (full.ans == 1 && dtype %in% c(5, 15)) 
        cat(" \nnote that the AIC of a fitted model for continuous data with litter effects \n  cannot be compared to the full model\n")
    if (output && model.sel != 1) {
        cat("\n ------------------------------------------------------- \n")
        cat("selected model:", model.txt[row.sel], " \n\n")
        if (aic.tmp != 0) 
            cat("using crit diff in AIC of", aic.tmp, " \n\n")
    }
    if (output && ans.all.tmp$trend) {
        cat(ans.all.tmp$report.pars)
        cat("------------------------------------------------------ \n\n")
    }
    ans.all.tmp <- f.move.sublist(ans.all, ans.all.spec)
    if (dtype %in% c(1, 5, 15) && !ans.all$MA.running) 
        ans.all.tmp <- f.resid.con(ans.all.tmp)
    if (dtype %in% c(5, 15)) {
        if (max(ans.all.tmp$Vdetlim, na.rm = T) != 0) 
            cat("\n litter effects not implemented for nonzero detection limit\n ")
        if (dtype == 5 && ans.all.tmp$nest.no == 0) 
            dtype <- 1
        else {
            ans.all.tmp <- f.nested.con(ans.all.tmp)
            ans.all.tmp$inter.var <- ans.all.tmp$inter.var
            ans.all.tmp$intra.var <- ans.all.tmp$intra.var
            df1 <- ans.all.tmp$df1
            df2 <- ans.all.tmp$df2
            nest.size <- ans.all.tmp$nest.size
            no.of.nests <- ans.all.tmp$no.of.nests
        }
    }
    ans.all.tmp$factor.name <- f.factorname(ans.all.tmp)
    ans.all.tmp$model.txt <- model.txt
    ans.all.tmp$Vloglik <- Vloglik
    ans.all.tmp$Vnpar <- Vnpar
    ans.all.tmp$Vaic <- Vaic
    ans.all.tmp$Vconverged <- Vconverged
    ans.all.tmp$MLE.all <- MLE.all
    ans.all.tmp$regr.par.matr <- f.pars(ans.all.tmp)$regr.par.matr
    ans.all.tmp$row.sel <- row.sel
    ans.all.tmp$conf.int <- matrix(rep(NA, 2), ncol = 2)
    if (cont) 
        ans.all.tmp$CED <- ans.all.tmp$regr.par[(ans.all.tmp$nr.aa + 
            1):(ans.all.tmp$nr.aa + ans.all.tmp$nr.bb)]
    if (model.sel != 1) {
        if (ans.all$dtype %in% c(5, 15)) 
            if (ans.all$do.MA) 
                ans.all.tmp$no.CI <- TRUE
            else ans.all.tmp$no.CI <- FALSE
        if (cont) 
            f.hit.constr(ans.all.tmp)
        if (!ans.all$no.CI && dtype == 3) 
            ans.all.tmp <- f.CI.sel.ord(ans.all.tmp)
        else if (ans.all$NES.ans == 2) {
            ans.all.tmp <- f.refit.nes(ans.all.tmp, interactive_mode = interactive_mode, results_env = results_env)
            # assign(ans.all.tmp$modelname, ans.all.tmp, envir = results_env)
            current_model <- ans.all.tmp$modelname
            ans.all$do.MA <- ans.all.tmp$do.MA
            ans.all$CES <- ans.all.tmp$CES
            # assign(current_model, ans.all.tmp, envir = results_env) # well, that's not right either
            # assign(current_model, ans.all, envir = results_env) 
            if (output) {
                cat(ans.all.tmp$report.pars)
                cat("------------------------------------------------------\n ")
            }
        }
        else if (!ans.all$no.CI && cont) {
            ans.all.tmp <- f.CI.sel(ans.all.tmp, results_env = results_env, interactive_mode = interactive_mode)
        }
        else if (!cont) {
            ans.all.tmp <- f.CI.sel(ans.all.tmp, results_env = results_env, interactive_mode = interactive_mode)
        }
    }
    if (model.sel == 1 && dtype %in% c(4, 6)) {
        ans.all.tmp$CED.matr <- matrix(rep(NA, 2), ncol = 2)
        ans.all.tmp$rank.low <- 1
    }
    if (model.switch == 0 && dtype %in% c(2, 4, 6)) 
        ans.all$ced.start.bin <- ans.all.tmp$CED
    loglik.summary <- data.frame(model = ans.all.tmp$model.txt, 
        converged = ans.all.tmp$Vconverged, loglik = ans.all.tmp$Vloglik, 
        npar = ans.all.tmp$Vnpar, AIC = ans.all.tmp$Vaic)
    ans.all.tmp$loglik.summary <- loglik.summary
    if (model.switch == 0) 
        ans.all[["EXP"]] <- f.model.specific.results(ans.all.tmp)
    if (model.switch == 1) 
        ans.all[["HILL"]] <- f.model.specific.results(ans.all.tmp)
    if (model.switch == 2) 
        ans.all[["INVEXP"]] <- f.model.specific.results(ans.all.tmp)
    if (model.switch == 3) 
        ans.all[["LOGN"]] <- f.model.specific.results(ans.all.tmp)
    if (model.switch == 0) 
        ans.all$logliks <- ans.all$EXP$loglik.summary
    if (model.switch == 1) 
        ans.all$logliks <- rbind(ans.all$EXP$loglik.summary, 
            ans.all$HILL$loglik.summary)
    if (model.switch == 2) 
        ans.all$logliks <- rbind(ans.all$EXP$loglik.summary, 
            ans.all$HILL$loglik.summary, ans.all$INVEXP$loglik.summary)
    if (model.switch == 3) 
        ans.all$logliks <- rbind(ans.all$EXP$loglik.summary, 
            ans.all$HILL$loglik.summary, ans.all$INVEXP$loglik.summary, 
            ans.all$LOGN$loglik.summary)
    ans.all$npar.null <- ans.all$NULL.model$npar
    ans.all$npar.aic.null <- ans.all$NULL.model$npar.aic
    ans.all$loglik.null <- ans.all$NULL.model$loglik
    ans.all$aic.null <- ans.all$NULL.model$aic
    if (length(ans.all$FULL.model) > 0) {
        ans.all$npar.full <- ans.all$FULL.model$npar
        ans.all$npar.aic.full <- ans.all$FULL.model$npar.aic
        ans.all$loglik.full <- ans.all$FULL.model$loglik
        ans.all$aic.full <- ans.all$FULL.model$aic
    }
    ans.all$fitted <- T
    if (model.switch == 0) {
        ans.all$row.sel.m1 <- row.sel.m1
        ans.all$aic.m1 <- aic.m1
        if (covar.no > 0) {
            ans.all$vv <- vv
        }
        if (cont) {
            ans.all$nr.var <- nr.var
            ans.all$var.start <- var.start
        }
        ans.all$Vaa.start <- Vaa.start
        ans.all$aa.start <- aa.start
        ans.all$Vced.start <- Vced.start
        ans.all$ced.start <- ced.start
        if (length(xans) > 1) 
            ans.all$RPF.start <- RPF.start
        ans.all$cc.start <- cc.start
        if (!cont) 
            ans.all$ths.start <- ths.start
        ans.all$skip.b <- skip.b
        ans.all$skip.0 <- skip.0
        ans.all$nr.first.models <- nr.first.models
        ans.all$model.type <- 2
        if (cont && !(model.sel == 1)) 
            if (max(ans.all$EXP$fct3) != ans.all$EXP$nr.var) {
                print("error in f.select.m5.con")
                print(ans.all$EXP$fct3)
                print(ans.all$EXP$nr.var)
                f.press.key.to.continue()
            }
    }
    if (0) 
        if (!(model.sel == 1 & model.switch == 1)) {
            if (max(fct1) > 1) 
                ans.all$fct1.no <- covar.no
            if (max(fct2) > 1) 
                ans.all$fct2.no <- covar.no
            if (max(fct3) > 1) 
                ans.all$fct3.no <- covar.no
        }
    ans.all$CED.model <- TRUE
    if (dtype %in% c(5, 15)) {
        ans.all$inter.var <- ans.all.tmp$inter.var
        ans.all$intra.var <- ans.all.tmp$intra.var
    }
    if (ans.all$dtype %in% c(5, 15)) {
        ans.all$df1 <- df1
        ans.all$df2 <- df2
        ans.all$no.of.nests <- no.of.nests
        ans.all$nest.size <- nest.size
    }
    #assign(current_model, ans.all, envir = results_env) # Nope... this only gets one model in each "loop" of the code somehow.
    if (exists("track")) 
        print("f.select.m5.con:   END")
    return(ans.all)
}



f.graph.window <- function(nr.pl = 1, nr.gr = 1, WAPP = FALSE, title = "", name.wapp = NA, 
    plotprefix = NA, quick.ans, svg.plots = FALSE, output = TRUE, filename = NULL, output_type = NULL) {
    if (exists("track")) 
        print("f.graph.window")
    if (nr.pl > 64) 
        cat("\nATTENTION:  nr of plots too large to plot together\n")
    if (!exists(".gw.size")) 
        .gw.size <- c(5, 5)
    if (WAPP && !svg.plots) 
        .gw.size <- c(400, 400)
    else if (WAPP && svg.plots) 
        .gw.size <- c(6, 6)
    f.create.graphwin(.gw.size[1], .gw.size[2], nr.gr = nr.gr, 
        WAPP = WAPP, title = title, name.wapp = name.wapp, plotprefix = plotprefix, 
        svg.plots = svg.plots, output = output, filename = filename, output_type = output_type)
    if (nr.pl == 1) {
        par(mfcol = c(1, 1))
        par(mex = 0.7)
        par(mar = c(6, 5, 6, 13))
    }
    if (nr.pl %in% 2) {
        par(mfcol = c(2, 1))
        par(mar = c(5.5, 4.8, 4, 12))
        par(cex.main = 1.2)
        par(cex.sub = 1)
        par(cex.lab = 1)
        par(cex = 0.6)
    }
    if (nr.pl %in% 3:4) {
        par(mfcol = c(2, 2))
        par(mar = c(3.2, 3.7, 2, 1.5))
        par(cex.main = 1.2)
        par(cex.sub = 0.9)
        par(cex.lab = 0.9)
        par(cex = 0.9)
    }
    if (nr.pl %in% 5:6) {
        par(mfcol = c(3, 2))
        par(mar = c(5, 3.5, 2.5, 2.5))
        par(cex.main = 1.2)
        par(cex.sub = 0.6)
        par(cex.lab = 1.1)
        par(cex = 0.6)
    }
    if (nr.pl %in% 7:9) {
        par(mfcol = c(3, 3))
        par(mar = c(2.5, 3, 2.5, 0.5))
        par(cex.main = 1.2)
        par(cex.sub = 0.6)
        par(cex.lab = 1.3)
        par(cex = 0.6)
    }
    if (nr.pl %in% 10:12) {
        par(mfcol = c(4, 3))
        par(mar = c(2.5, 2.5, 1.5, 1.5))
        par(cex.main = 1.2)
        par(cex.sub = 0.6)
        par(cex.lab = 0.6)
        par(cex = 0.6)
    }
    if (nr.pl %in% 13:16) {
        par(mfcol = c(4, 4))
        par(mar = c(2, 2, 1.2, 0.5))
    }
    if (nr.pl %in% 17:20) {
        par(mfcol = c(4, 5))
        par(mar = c(3.3, 2.5, 2, 1.5))
    }
    if (nr.pl %in% 21:25) {
        par(mfcol = c(5, 5))
        par(mar = c(3, 2, 1.8, 1.2))
    }
    if (nr.pl %in% 26:30) {
        par(mfcol = c(6, 5))
        par(mar = c(2, 1.5, 1.5, 1))
    }
    if (nr.pl %in% 31:36) {
        par(mfcol = c(6, 6))
        par(mar = c(2, 1.5, 1.5, 1))
    }
    if (nr.pl %in% 37:42) {
        par(mfcol = c(7, 6))
        par(mar = c(2, 1.5, 1.5, 1))
    }
    if (nr.pl %in% 43:49) {
        par(mfcol = c(7, 7))
        par(mar = c(2, 1.5, 1.5, 1))
    }
    if (nr.pl %in% 50:56) {
        par(mfcol = c(8, 7))
        par(mar = c(2, 1.5, 1.5, 1))
    }
    if (nr.pl %in% 57:64) {
        par(mfcol = c(8, 8))
        par(mar = c(2, 1.5, 1.5, 1))
    }
    if (exists("track")) 
        print("f.graph.window   END")
    return(invisible())
}



f.create.graphwin <- function(aa, bb, title = "", name.wapp = NA, WAPP = FALSE, plotprefix = NA, 
    nr.gr = 1, svg.plots = FALSE, output = TRUE, output_type = NULL, filename = NULL) {
    if (exists("track")) 
        print("f.create.graphwin")
    lst <- dev.list()
    if (length(lst) > 60) {
        cat("\n ATTENTION:  too many graphical windows are opened")
        cat("\nthey will be deleted, this is the last chance to save any of them")
        if (output) 
            f.press.key.to.continue("\npress any key if they can be deleted ... > ")
        f.delete.gw()
        assign(".ypos", 0)
    }
    if (nr.gr > 25) 
        f.delete.gw()
    if (WAPP) {
        if (!is.null(filename)) {
            filename <- paste0(filename, "_", name.wapp)
        }
        if (svg.plots & is.null(filename)) {
            svg(paste(plotprefix, name.wapp, ".svg", sep = ""), 
                aa, bb)
        } else if (output_type == "jpeg") {
            filename <- paste0(filename, ".jpeg")
            grDevices::jpeg(filename, width = aa, height = bb)
        } else if (output_type == "pdf") {
            filename <- paste0(filename, ".pdf")
            grDevices::pdf(filename, width = aa, height = bb)
        } else if (output_type == "png") {
            filename <- paste0(filename, ".png")
            grDevices::png(filename, width = aa*3, height = bb*3, res = 300)
        } else if (output_type == "tiff") {
            filename <- paste0(filename, ".tiff")
            grDevices::tiff(filename, width = aa*3, height = bb*3, res = 300)
        } else if (output_type == "svg") {
            message("Making SVG...")
            message(filename)
            filename <- paste0(filename, ".svg")
            grDevices::svg(filename, width = aa, height = bb)
        } else if (output_type == "none") {
            message("Not saving PROAST plots to file!")
        } else {
            png(paste(plotprefix, name.wapp, ".png", sep = ""), 
            aa, bb)
        }
    }
    else {
        title <- paste(dev.cur(), "PROAST:", title)
        y.pos <- .ypos + 5
        if (y.pos > 100) 
            y.pos <- 0
        assign(".ypos", y.pos, immediate = T, pos = 1)
        if (.Platform$OS.type == "windows") 
            windows(aa, bb, 10, ypos = y.pos, restoreConsole = T, 
                title = title)
        else X11(, aa, bb, 10, title = title, ypos = y.pos)
    }
    if (exists("track")) 
        print("f.create.graphwin:   END")
    return(invisible())
}



f.qfit <- function(ans.all.fit, plot.type, plt.mns = 3, output = TRUE, display_plots = TRUE) {
    if (exists("track")) 
        print("f.qfit")
    plot.type.save <- ans.all.fit$plot.type
    ans.all.fit$plot.type <- plot.type
    ans.all.fit$plt.mns <- plt.mns
    if (!output) 
        ans.all.fit$plot.type <- 0
    with(ans.all.fit, {
        if (cont) {
            increase.00 <- ans.all.fit$increase
            ans.all.fit$text.par <- f.text.par(ans.all.fit)
            par.start.tmp <- ans.all.fit$par.start
            if (any(!is.finite(par.start.tmp))) {
                ans.all.fit <- f.start.con(ans.all.fit, adjust = F, 
                  fitted = F, tmp.quick = T)
                if (ans.all.fit$adjust.start && !WAPP) {
                  cat("\nadjust start values")
                  ans.all.fit <- f.start.con(ans.all.fit, adjust = T, 
                    fitted = F, tmp.quick = T)
                  par.start <- ans.all.fit$par.start
                }
                if (ans.all.fit$par.start[1] == 0) {
                  ans.all.fit <- f.start.con(ans.all.fit, adjust = T)
                }
            }
            else {
                ans.all.fit$regr.par <- par.start.tmp[-(1:max(fct3))]
                ans.all.fit$par.start <- par.start.tmp
                ans.all.fit$npar <- length(par.start.tmp)
                ans.all.fit <- f.constr.con(ans.all.fit, tmp.quick = T)
            }
            ans.all.fit$increase <- increase.00
            ans.all.fit <- f.nlminb(ans.all.fit, tmp.quick = T)
            if (ans.all.fit$MLE[1] == 0) 
                return(ans.all.fit)
            ans.all.fit <- f.pars(ans.all.fit)
            if (output) 
                ans.all.fit$report.pars <- f.report.pars(ans.all.fit)
            MLE <- ans.all.fit$MLE
            ans.all.fit$regr.par <- MLE[(max(fct3) + 1):length(MLE)]
            if (length(xans) > 1) {
                RPF.vec <- ans.all.fit$regr.par[(length(ans.all.fit$regr.par) - 
                  nr.dosecol + 2):length(ans.all.fit$regr.par)]
                if (model.ans != 42) 
                  x <- Mx %*% c(1, RPF.vec)
                ans.all.fit$x <- x
                ans.all.fit$xy.lim[2:3] <- c(min(x), max(x))
                if (min(x) == 0) 
                  dum.contr <- min(x[x > 0])/4
                else dum.contr <- min(x)
                ans.all.fit$xy.lim[1] <- dum.contr
            }
            if (interrupt && plot.type > 0 && !WAPP) {
                ans.all.fit$heading <- ""
                ans.all.fit$xy.lim[2:3] <- c(min(x), max(x))
                ignore <- f.plot.con(ans.all.fit, display_plots = display_plots)
                modelname <- ans.all.fit$modelname
                if (display_plots == TRUE) {
                  title(main = paste(modelname))
                }
                ans.all.fit$l.ty = 1
                if (model.ans != 11) 
                  f.lines.con(ans.all.fit, display_plots = display_plots)
            }
        }
        if (!cont) {
            ans.all.fit$nr.aa <- max(fct1)
            ans.all.fit$nr.bb <- max(fct2)
            if (model.type == 1) {
                if (model.ans != 14) {
                  cat("\n error in f.qfit: model.type and model.ans do not match\n")
                  stop()
                }
                ans.all.fit <- f.start.bin(ans.all.fit, tmp.quick = T)
                ans.all.fit <- f.nlminb(ans.all.fit, tmp.quick = T)
                if (dtype == 6) {
                  ans.all.fit$pi.full <- ans.all.fit$MLE[-1]
                  ans.all.fit$alfa.start <- ans.all.fit$MLE[1]
                }
                if (interrupt && plot.type > 0 && hill == 0) {
                  ans.all.fit$heading <- ""
                  ans.all.fit$xy.lim[2:3] <- c(min(x), max(x))
                  ans.all.fit$shift.tmp <- f.plot.frq(ans.all.fit)$shift.tmp
                  title(main = paste(modelname))
                }
                ans.all.fit$npar <- length(ans.all.fit$MLE)
            }
            if (model.type == 2) {
                par.start.tmp <- ans.all.fit$par.start
                if (any(is.na(par.start.tmp)) || length(xans) > 
                  1) {
                  ans.all.fit <- f.start.cat(ans.all.fit, tmp.quick = T)
                  nrp <- ans.all.fit$nrp
                  ans.all.fit$text.par <- f.text.par(ans.all.fit)
                }
                else {
                  nrp <- length(par.start.tmp) - nth - 1
                  if (dtype == 6) 
                    nrp <- nrp - 1
                  loglik.start <- -f.lik.cat(par.start.tmp, x, 
                    y, kk, nn, dtype, fct1, fct2, nrp, nth, nr.aa, 
                    nr.bb, model.ans, model.type, CES.cat = CES.cat, 
                    ttt = ttt, twice = twice, CES = CES, decr.zz = decr.zz, 
                    cc.inf = cc.inf, cens.up = cens.up, fct3 = fct3, 
                    ces.ans = ces.ans, Mx = Mx, kk.tot = kk.tot, 
                    nn.tot = nn.tot, quick.ans = quick.ans)
                  if (is.na(loglik.start) | loglik.start == -1e+12) {
                    ans.all.fit <- f.start.cat(ans.all.fit, tmp.quick = T)
                  }
                  else {
                    ans.all.fit$par.start <- par.start.tmp
                    ans.all.fit$npar <- length(par.start.tmp)
                  }
                  ans.all.fit$text.par <- f.text.par(ans.all.fit)
                }
                ans.all.fit$nrp <- nrp
                ans.all.fit <- f.constr.con(ans.all.fit, tmp.quick = T)
                ans.all.fit <- f.nlminb(ans.all.fit, tmp.quick = T)
                par.lst <- f.split.par(ans.all.fit$MLE, ans.all.fit$nrp, 
                  ans.all.fit$nth, dtype, fct3)
                regr.par <- par.lst$regr.par
                th.par <- par.lst$th.par
                sig.par <- par.lst$sig.par
                ans.all.fit <- f.pars(ans.all.fit)
                ans.all.fit$report.pars <- f.report.pars(ans.all.fit)
                if (interrupt && plot.type > 0 && hill == 0) {
                  print("f.qfit 000000")
                  print("Attention to programmer: this part of the code is used, which is unexpected")
                  f.press.key.to.continue()
                  ans.all.fit$heading <- ""
                  ans.all.fit$xy.lim[2:3] <- c(min(x), max(x))
                  ans.all.fit$regr.par <- regr.par
                  ans.all.fit$th.par <- th.par
                  ans.all.fit$sig.par <- sig.par
                  if (dtype %in% 2:3) {
                    ans.all.fit$plot.type <- plot.type.save
                    f.lines.cat(ans.all.fit)
                    ans.all.fit$regr.par <- regr.par
                  }
                  else {
                    kk.dum <- 1
                    f.plot.frq(ans.all.fit)
                    f.lines.frq(ans.all.fit)
                    ans.all.fit$regr.par <- regr.par
                  }
                  title(main = paste(modelname))
                }
            }
        }
        if (!is.finite(ans.all.fit$loglik)) 
            ans.all.fit$loglik <- -1e+12
        ans.all.fit$plot.type <- plot.type.save
        ans.all.fit$fitted <- TRUE
        if (exists("track")) 
            print("f.qfit:  END")
        return(ans.all.fit)
    })
}



f.text.par <- function(ans.all, brief = F) {
    if (exists("track")) 
        print("f.text.par")
    with(ans.all, {
        if (brief) {
            nr.var <- 1
            nr.aa <- 1
            nr.bb <- 1
            nr.cc <- 1
            nr.dd <- 1
            nr.RPF <- 1
            if (incr.decr.no > 0) 
                nr.cc <- 2
        }
        if (!cont) 
            nr.var <- 0
        if (dtype == 6) 
            nr.var <- max(fct3)
        if (nr.aa == 1) {
            fct1.txt <- NULL
        }
        CED.name <- "CED"
        if (dtype %in% c(4, 6) && ces.ans > 1) 
            CED.name <- "BMD"
        if (nr.bb == 1 && length(xans) == 1 && model.ans != 46) {
            fct2.txt <- NULL
        }
        if (nr.var == 1) {
            fct3.txt <- NULL
        }
        if (nr.cc == 1) {
            fct4.txt <- NULL
        }
        if (nr.dd == 1) {
            fct5.txt <- NULL
        }
        if (!cont && model.type == 1) {
            switch(model.ans, text.par <- paste("a", fct1.txt[1:nr.aa], 
                sep = "-"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c"), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c", 
                  "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c"), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c", 
                  "dd"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), text.par <- c(rep("aa", nr.aa), rep("bb", 
                  nr.bb)), text.par <- c(rep("aa", nr.aa), rep("bb", 
                  nr.bb), "cc"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("bb", fct2.txt[1:nr.bb], 
                  sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
                {
                  if (length(ans.all$x.full) == 0) x.full <- x
                  if (dtype %in% c(4, 6, 84)) text.par <- paste("group", 
                    1:length(x.full))
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-")) else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"))
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"), "c") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"), "c")
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"), "c", "d") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"), "c", "d")
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"), "c") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"), "c")
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"), "c") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"), "c")
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"), "c", "d") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"), "c", "d")
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"), "c") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"), "c")
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"), "c") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"), "c")
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-")) else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"))
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-"), "c") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"), "c")
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-")) else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"))
                }, {
                  if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                    fct2.txt[1:nr.bb], sep = "-")) else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                    fct2.txt[1:nr.bb], sep = "-"))
                }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "BMDratio"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"))
            if (!(model.ans %in% c(14))) 
                if (length(xans) > 1) {
                  text.RPF <- character(0)
                  for (ii in 2:nr.dosecol) text.RPF <- c(text.RPF, 
                    paste("RPF", fct2.txt[ii], sep = "-"))
                  text.par <- c(text.par, text.RPF)
                }
        }
        if (cont || model.type == 2) {
            switch(model.ans, text.par <- paste("a", fct1.txt[1:nr.aa], 
                sep = "-"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                    fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-")), {
                  if (ans.m6.sd == 1) {
                    text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-"), paste("c", fct4.txt[1:nr.cc], 
                      sep = "-"), paste("d", fct5.txt[1:nr.dd], 
                      sep = "-"))
                  }
                  if (ans.m6.sd == 2) text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), "b", "q", 
                    "d ")
                }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                    fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-")), text.par <- rep("GM", 
                  length(regr.par)), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste(CED.name, 
                  fct2.txt[1:nr.bb], sep = "-")), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste(CED.name, 
                  fct2.txt[1:nr.bb], sep = "-"), paste("d", fct5.txt[1:nr.dd], 
                  sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste(CED.name, fct2.txt[1:nr.bb], 
                  sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste(CED.name, fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), "c", "BMDratio"), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("d", fct5.txt[1:nr.dd], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                    fct4.txt[1:nr.cc], sep = "-")), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-"), 
                  paste("d", fct5.txt[1:nr.dd], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                    fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste(CED.name, 
                  fct2.txt[1:nr.bb], sep = "-")), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste(CED.name, 
                  fct2.txt[1:nr.bb], sep = "-"), paste("d", fct5.txt[1:nr.dd], 
                  sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste(CED.name, fct2.txt[1:nr.bb], 
                  sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste(CED.name, fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), "c"), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), "c"), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), "c", "d"), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), "c", "d"), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), "c"), {
                  if (brief) text.par <- c("a", "b1", "b2", "c1", 
                    "c2", "dd") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), "bb1", paste("bb2", 
                    fct2.txt[1:nr.bb], sep = "-"), paste("cc1", 
                    fct4.txt[1:nr.cc], sep = "-"), "cc2", "dd")
                }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("bb1", fct2.txt[1:nr.bb], 
                  sep = "-"), "bb2", "cc", "dd1", "dd2"), {
                  if (brief) text.par <- c("a", "b1", "b2", "c1", 
                    "c2", "d1", "d2") else text.par <- c(paste("a", 
                    fct1.txt[1:nr.aa], sep = "-"), "bb1", paste("bb2", 
                    fct2.txt[1:nr.bb], sep = "-"), paste("cc1", 
                    fct4.txt[1:nr.cc], sep = "-"), "cc2", "d1", 
                    "d2")
                }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "cc", "dd"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c", 
                  "d"), , , , text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  "b1", "b2", "c", "d"), text.par <- c(paste("a", 
                  fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                  sep = "-"), "c", paste("d", fct5.txt[1:nr.dd], 
                  sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c1", 
                  "c2", "d"), {
                  if (brief) text.par <- c("a", "CED/RPF", "c", 
                    "d") else {
                    nr.bb <- length(fct2.txt)
                    ced.txt <- paste("RPF", fct2.txt[1:nr.bb], 
                      sep = "-")
                    ced.txt[ref.lev] <- paste(CED.name, fct2.txt[ref.lev], 
                      sep = "-")
                    text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), ced.txt, paste("c", fct4.txt[1:nr.cc]), 
                      paste("d", fct5.txt[1:nr.dd], sep = "-"))
                  }
                }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), "GCED", "q", "d "))
            if (model.ans == 48 || (dtype == 4 && model.type == 
                2 && model.ans == 34)) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"), "aa.t", "bb.t", 
                  "dd.t")
            if (model.ans == 50) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"), "aa.t", "bb.t", 
                  "dd.t")
            if (model.ans == 49) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("ED50", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans == 51) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste(CED.name, fct2.txt[1:nr.bb], sep = "-"), 
                  paste("d", fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans == 52) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste(CED.name, fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans == 53) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste(CED.name, fct2.txt[1:nr.bb], sep = "-"), 
                  paste("d", fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans == 54) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste(CED.name, fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans %in% 55:56) 
                text.par <- c(paste("Km", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("Vmax", fct2.txt[1:nr.bb], 
                  sep = "-"), paste("V", fct4.txt[1:nr.cc], sep = "-"))
            if (model.ans == 57) {
                if (brief) 
                  text.par <- c("a", "b1", "b2", "c1", "d1", 
                    "d2")
                else text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("bb1", fct5.txt[1:nr.dd], 
                  sep = "-"), paste("bb2", fct2.txt[1:nr.bb], 
                  sep = "-"), paste("cc1", fct4.txt[1:nr.cc], 
                  sep = "-"), "d1", "d2")
            }
            if (model.ans == 58) 
                text.par <- c("FD.V", "ka", "ke")
            if (model.ans == 59) 
                text.par <- c("C1", "C2", "ka", "lab1", "lab2")
            if (model.ans == 60) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("CED", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans == 61) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("CED", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans == 62) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                    fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans == 63) 
                text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("CED", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                    fct5.txt[1:nr.dd], sep = "-"))
            if (model.ans == 64) 
                text.par <- c("D", paste("V", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("lab.min", fct2.txt[1:nr.bb], 
                  sep = "-"), paste("lab.plus", fct4.txt[1:nr.cc], 
                  sep = "-"), paste("kc", fct5.txt[1:nr.dd], 
                  sep = "-"))
            if (model.ans == 65) 
                text.par <- c(paste("cc", fct1.txt[1:nr.aa], 
                  sep = "-"), paste("k1", fct2.txt[1:nr.bb], 
                  sep = "-"), paste("ke", fct5.txt[1:nr.dd], 
                  sep = "-"))
            if (cont && fit.ans == 1) 
                text.par <- c(paste("var", fct3.txt[1:nr.var], 
                  sep = "-"), text.par)
            if (!(model.ans %in% c(11, 42))) 
                if (length(xans) > 1) {
                  text.RPF <- character(0)
                  if (brief) 
                    text.par <- c(text.par, "RPF")
                  else {
                    for (ii in 2:nr.dosecol) text.RPF <- c(text.RPF, 
                      paste("RPF", fct2.txt[ii], sep = "-"))
                    text.par <- c(text.par, text.RPF)
                  }
                }
            if (!cont && ces.ans != 4) {
                if (model.ans %in% c(12:15, 22:25) && nr.aa > 
                  1 && nr.bb == 1) 
                  for (ii in 2:(nr.aa + 1)) {
                    text.par[ii] <- paste(CED.name, fct1.txt[ii - 
                      1], sep = "-")
                  }
            }
            if (!cont) {
                if (dtype == 4) 
                  text.th <- "th"
                else text.th <- paste("th-", 1:nth, sep = "")
                if (dtype == 3) 
                  text.par <- c(text.par, text.th, "sigma")
                if (dtype == 4) {
                  if (length(unique(fct3)) == 1) 
                    text.par <- c(text.par, text.th, "sigma")
                  if (max(fct3) > 1) 
                    text.par <- c(text.par, "sigma")
                }
            }
        }
        if (dtype == 6) 
            text.par <- c(paste("alfa", fct3.txt[1:nr.var], sep = "-"), 
                text.par)
        if (exists("track")) 
            print("f.text.par:  END")
        return(text.par)
    })
}



f.start.con <- function(ans.all, adjust = F, fitted = F, tmp.quick = F, display_plots = TRUE) {
    if (exists("track")) 
        print("f.start.con")
    if (ans.all$model.ans == 0) {
        cat("\n\nchoose model first\n")
        f.press.key.to.continue()
        return(ans.all)
    }
    if (ans.all$cont && ans.all$model.ans == 10 && ans.all$list.logic) {
        ans.dum <- menu(c("yes", "no"), title = "\ndo you want to use previous fitted model 5 for start values\n?")
        if (ans.dum == 1) {
            model.cur <- as.numeric(readline(prompt = "\ngive model number of previous fit >"))
            MLE <- ans.all$MLE
            nr.var <- ans.all$nr.var
            nr.aa <- ans.all$nr.aa
            nr.bb <- ans.all$nr.bb
            par.start <- MLE
            cc <- MLE[length(MLE) - 1]
            dd <- MLE[length(MLE)]
            tmp.bb <- MLE[(nr.var + nr.aa + 1):(nr.var + nr.aa + 
                nr.bb)]
            if (model.cur == 15) {
                tmp.bb <- f.bb.con(15, cc, dd, tmp.bb, ans.all$CES)
            }
            bb.new <- (1/tmp.bb)^(1/dd)
            par.start[(nr.var + nr.aa + 1):(nr.var + nr.aa + 
                nr.bb)] <- bb.new
            ans.all$par.start <- par.start
            ans.all$regr.par <- ans.all$par.start[-(1:nr.var)]
            ans.all <- f.plot.con(ans.all, display_plots = display_plots)
            f.lines.con(ans.all)
            with(ans.all, {
                loglik.start <- -f.lik.con(par.start, x, y, dtype, 
                  fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                  nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                  ttt = ttt, trace = F, fct4 = fct4, fct5 = fct5, 
                  cens.up = cens.up, par.tmp = NA, increase = increase, 
                  ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd, 
                  Mx = Mx, x1 = 1, x2 = x2, cc.inf = cc.inf)
                cat("\n\nlog-likelihood at start values: ", round(loglik.start, 
                  2), "\n\n")
            })
            ans.all <- f.constr.con(ans.all, tmp.quick = T)
            return(ans.all)
        }
    }
    if (!adjust) 
        if (ans.all$model.ans != 16 && !tmp.quick && max(ans.all$fct5) == 
            1) 
            ans.all <- f.clear(ans.all)
    ans.all$adjust <- adjust
    ans.all$fitted <- fitted
    ans.all$tmp.quick <- tmp.quick
    ans.all$nr.aa <- max(ans.all$fct1)
    ans.all$nr.bb <- max(ans.all$fct2)
    if (ans.all$model.ans == 42) 
        ans.all$nr.bb <- 2
    ans.all$nr.var <- max(as.numeric(ans.all$fct3))
    ans.all$nr.cc <- max(ans.all$fct4)
    ans.all$nr.dd <- max(ans.all$fct5)
    ans.all$heading <- "start values"
    with(ans.all, {
        if (fit.ans == 2) 
            nr.var <- 0
        cc.upp <- 10000
        dd.low <- 0.1
        dd.upp <- 10
        if (dtype == 5) 
            dtype <- 1
        trace.tmp <- F
        if (quick.ans == 1) 
            trace.tmp <- T
        y.sort <- yy[order(x)]
        x.sort <- sort(x)
        lower.var <- 1e-06
        upper.var <- 10
        if (dtype %in% c(25, 250, 26, 260)) 
            upper.var <- Inf
        if (adjust == F) {
            if (dtype %in% c(10, 15, 250, 260) & length(x) < 
                10) {
                xtmp <- x.sort
                ytmp <- y.sort
            }
            else {
                sub <- round(length(x.sort)/5)
                if (sub == 0) 
                  sub <- 1
                xtmp <- numeric()
                ytmp <- numeric()
                regr.par <- numeric()
                e1 <- 0
                e2 <- 0
                for (k in 1:5) {
                  e2 <- e2 + sub
                  qq <- x.sort[e1:e2]
                  xtmp[k] <- mean(qq[is.na(qq) == F])
                  qq <- y.sort[e1:e2]
                  if (dtype %in% c(25, 250)) 
                    ytmp[k] <- mean(qq[is.na(qq) == F])
                  if (dtype %in% c(26, 260)) 
                    ytmp[k] <- mean(sqrt(qq[is.na(qq) == F]))^2
                  else if (!(dtype %in% c(25, 26, 250, 260))) 
                    ytmp[k] <- exp(mean(logb(qq[is.na(qq) == 
                      F])))
                  e1 <- e1 + sub
                }
            }
            if (length(ytmp[!is.na(ytmp)]) == 3) {
                ytmp[4:5] <- ytmp[3]
                xtmp[4:5] <- xtmp[3]
            }
            if (length(ytmp[!is.na(ytmp)]) == 4) {
                ytmp[5] <- ytmp[4]
                xtmp[5] <- xtmp[4]
            }
            top <- length(ytmp)
            aa <- NA
            bb <- NA
            cc <- NA
            dd <- NA
            if (model.ans != 1) {
                lm.res <- f.start.lm.con(xtmp, ytmp, model.ans, 
                  dtype)
                aa <- abs(lm.res[1])
                if (model.ans %in% c(2, 3, 7, 8, 17, 18)) 
                  bb <- lm.res[2]
                else bb <- abs(lm.res[2])
                if (quick.ans == 1 && model.type == 1) {
                  increase <- lm.res[3]
                  if (model.ans %in% c(21)) 
                    increase <- -increase
                }
                if ((model.ans == 6 && max(fct4) > 1) || incr.decr.no > 
                  0) 
                  increase <- 0
                if (!WAPP) {
                  if (cont && model.ans %in% c(12:15, 22:25, 
                    46, 51:54) && quick.ans == 1 && increase == 
                    -1 && CES > 0) {
                    incr.ans <- menu(c("yes", "no"), title = "\n you entered a positive value for CES, \n                             is the function increasing?\n")
                    if (incr.ans == 1) 
                      increase <- 1
                    if (incr.ans == 2) 
                      CES <- -CES
                  }
                  if (cont && model.ans %in% c(12:15, 22:25, 
                    46, 51:54) && quick.ans == 1 && increase == 
                    1 && CES < 0) {
                    ces.ans <- menu(c("yes", "no"), title = "\n you entered a negative value for CES, is the function decreasing?\n")
                    if (ces.ans == 1) 
                      increase <- -1
                    if (ces.ans == 2) 
                      CES <- -CES
                  }
                }
                if (incr.decr.no > 0 && model.ans %in% c(0, 15, 
                  25, 52, 54, 46)) {
                  CES.abs <- abs(CES[1])
                  CES <- rep(CES.abs, length(x))
                  CES[fct4 == 1] <- -CES.abs
                  ans.all$CES <- CES
                }
                if (incr.decr.no == 0) {
                  if (CES <= -1) {
                    CES.old <- CES
                    FC <- abs(CES) + 1
                    CES <- 1/FC - 1
                    if (output) {
                      cat("ATTENTION: you entered a value of", 
                        CES.old, "for CES, this is translated into", 
                        CES, "\n")
                      f.press.key.to.continue()
                    }
                  }
                }
                else {
                  cc <- c(1e-18, 1e+18)
                  nr.cc <- 1
                }
                ans.all$CES <- CES
            }
            if (incr.decr.no > 0) {
                cc <- c(1e-18, 1e+18)
                nr.cc <- 1
            }
            switch(model.ans, {
                if (nr.bb > 1) {
                  cat("\nThis model is not defined for different b parameters!")
                  f.press.key.to.continue()
                  ans.all$fct2 <- rep(1, length(x))
                }
                if (nr.cc > 1) {
                  cat("\nThis model is not defined for different c parameters!")
                  return(ans.all)
                }
                if (nr.dd > 1) {
                  cat("\nThis model is not defined for different d parameters!")
                  return(ans.all)
                }
                aa <- mean(ytmp)
                regr.par <- rep(aa, nr.aa)
            }, {
                if (fitted == T) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                }
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
            }, {
                if (fitted) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                }
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(dd, nr.dd))
            }, {
                if (fitted) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  cc <- ytmp[top]/aa
                } else {
                  cc <- ytmp[top]/aa
                  if (!is.finite(bb)) bb <- xtmp[5]/100
                }
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc))
                if (cc < 1) increase <- -1 else increase <- 1
            }, {
                if (fitted) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  cc <- ytmp[top]/aa
                } else {
                  cc <- ytmp[top]/aa
                  if (!is.finite(bb)) bb <- xtmp[5]/100
                  dd <- 1
                }
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
                if (cc < 1) increase <- -1 else increase <- 1
            }, {
                bb <- 1/bb
                dd <- 1
                if (fitted) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  cc <- ytmp[top]/aa
                  if (!cont) bb <- bb * 10
                } else {
                  cc <- ytmp[top]/aa
                  cc <- 1.01 * ytmp[top]/aa
                }
                if (!is.finite(bb)) bb <- xtmp[5]/100
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                if (fitted) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                } else bb <- 1/bb
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
            }, {
                if (fitted) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                } else bb <- 1/bb
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(dd, nr.dd))
            }, {
                if (fitted) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  cc <- ytmp[top]/aa
                } else {
                  bb <- abs(1/bb)
                  cc <- ytmp[top]/aa
                }
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc))
            }, {
                cc <- ytmp[top]/aa
                bb <- 1/bb
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                regr.par <- f.regr.par.full(ans.all)
                ans.all$regr.par <- regr.par
            }, {
                if (fitted == T) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                }
                if (fitted == F || is.na(bb) || bb == 0) {
                  if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                  CED <- (1/bb) * logb(CES + 1)
                  CED <- abs(CED)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                    min(x))/3
                }
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb))
            }, {
                if (fitted == T) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                }
                if (fitted == F || is.na(bb) || bb == 0) {
                  if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                  CED <- (1/bb) * logb(CES + 1)
                  CED <- abs(CED)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                    min(x))/3
                }
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(dd, nr.dd))
            }, {
                if (fitted == T) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  cc <- 0
                }
                if (fitted == F || is.na(bb) || bb == 0) {
                  cc <- ytmp[top]/aa
                  if (cont) {
                    if (CES < 0 & cc > 1 + CES) {
                      cc <- 1 + CES - (1 + CES)/100
                    }
                    if (CES > 0 & cc < 1 + CES) {
                      cc <- 1 + CES + CES/100
                    }
                  }
                  if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                  if (cc == 0) cc <- 0.1
                  if (cc == 1) cc <- 0.9
                  CED <- -(1/bb) * logb((CES + 1 - cc)/(1 - cc))
                  CED <- abs(CED)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                    min(x))/3
                }
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc))
            }, {
                dd <- 1
                if (incr.decr.no > 0) CED <- mean(xtmp) else {
                  cc <- ytmp[top]/aa
                  if (cont) {
                    if ((CES < 0) && (cc > 1 + CES)) {
                      cc <- 1 + CES - (1 + CES)/100
                    }
                    if ((CES > 0) & (cc < 1 + CES)) {
                      cc <- 1 + CES + CES/100
                    }
                    if (cc == 0) cc <- 0.1
                    if (cc == 1) cc <- 0.9
                    if (!cont) cc <- 0.1
                  }
                  CED <- (-(1/bb) * logb((-CES + 1 - cc)/(1 - 
                    cc)))^(1/dd)
                  CED <- abs(CED)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                    min(x))/3
                  if (is.na(CED) || CED == 0) CED <- mean(xtmp)
                }
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                regr.par <- MLE[-(1:nr.var)]
                if (!exists("CES.16")) {
                  cat("\nCES.16 does not exist and is created (0.10 and 0.05)\n")
                  CES.16 <- list(CES2 = 0.1, CES1 = 0.05)
                  f.assign("CES.16", CES.16)
                }
                CES2 <- CES.16$CES2
                CES1 <- CES.16$CES1
                aa.vec <- regr.par[1:nr.aa]
                bb.vec <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
                cc <- regr.par[nr.aa + nr.bb + 1]
                dd <- regr.par[length(regr.par)]
                if (cc < 1) {
                  CES2 <- -CES2
                  CES1 <- -CES1
                  if (cc > 1 + CES2) {
                    cat("\nATTENTION: value of c does not allow BMD ratio")
                    f.press.key.to.continue()
                    return(ans.all)
                  }
                }
                if (cc > 1) {
                  if (cc < 1 + CES2) {
                    cat("\nATTENTION: value of c does not allow BMD ratio")
                    return(ans.all)
                  }
                }
                BMDratio <- (log((CES2 + 1 - cc)/(1 - cc))/log((CES1 + 
                  1 - cc)/(1 - cc)))^(1/dd)
                regr.par <- c(aa.vec, bb.vec, rep(cc, nr.cc), 
                  BMDratio)
            }, {
                if (bb == Inf) bb <- 1000
                if (model.type == 2) increase <- -1
                if (increase == 1) bb <- -abs(bb)
                if (increase == -1) bb <- abs(bb)
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
                if (cont) {
                  y.expect <- f.expect.con(model.ans, x, regr.par, 
                    CES = CES, increase = increase, fct1 = fct1, 
                    fct2 = fct2, cc.inf = cc.inf)
                  y.expect <- y.expect[!is.na(y.expect)]
                  count <- 0
                  while (any(y.expect < 0) && count < 20) {
                    bb <- 1.5 * bb
                    regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                      rep(dd, nr.dd))
                    y.expect <- f.expect.con(model.ans, x, regr.par, 
                      CES = CES, increase = increase, fct1 = fct1, 
                      fct2 = fct2, cc.inf = cc.inf)
                    count <- count + 1
                  }
                }
            }, {
                if (bb == Inf) bb <- 1000
                if (model.type == 2) increase <- -1
                if (quick.ans == 1 && cont) if (bb < 0) increase <- 1 else increase <- -1
                if (increase == 1) bb <- -abs(bb)
                if (increase == -1) bb <- abs(bb)
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(dd, nr.dd))
                if (cont) {
                  y.expect <- f.expect.con(model.ans, x, regr.par, 
                    CES = CES, increase = increase, fct1 = fct1, 
                    fct2 = fct2, cc.inf = cc.inf)
                  y.expect <- y.expect[!is.na(y.expect)]
                  count <- 0
                  while (any(y.expect < 0) && count < 20) {
                    bb <- 1.5 * bb
                    regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                      rep(dd, nr.dd))
                    y.expect <- f.expect.con(model.ans, x, regr.par, 
                      CES = CES, increase = increase, fct1 = fct1, 
                      fct2 = fct2, cc.inf = cc.inf)
                    count <- count + 1
                  }
                }
            }, {
                cc <- ytmp[top]/aa
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc))
            }, {
                cc <- ytmp[top]/aa
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                cc <- ytmp[top]/aa
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                if (CES > 0) bb <- -bb
                CED <- f.inv.con(model.ans = 17, c(aa, bb), CES)
                if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                  min(x))/3
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb))
                if (cont) {
                  y.expect <- f.expect.con(model.ans, x, regr.par, 
                    CES = CES, increase = increase, fct1 = fct1, 
                    fct2 = fct2, cc.inf = cc.inf)
                  count <- 0
                  while (any(y.expect < 0) && count < 20) {
                    CED <- 1.5 * CED
                    regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb))
                    y.expect <- f.expect.con(model.ans, x, regr.par, 
                      CES = CES, increase = increase, fct1 = fct1, 
                      fct2 = fct2, cc.inf = cc.inf)
                    count <- count + 1
                  }
                }
            }, {
                if (CES > 0) bb <- -bb
                if (CES < 0) increase <- -1 else increase <- 1
                dd <- 1
                CED <- f.inv.con(model.ans = 18, c(aa, bb, dd), 
                  CES)
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(dd, nr.dd))
                if (cont) {
                  y.expect <- f.expect.con(model.ans, x, regr.par, 
                    CES = CES, increase = increase, fct1 = fct1, 
                    fct2 = fct2, cc.inf = cc.inf)
                  count <- 0
                  while (any(y.expect < 0) && count < 20) {
                    CED <- 1.5 * CED
                    regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                      rep(dd, nr.dd))
                    y.expect <- f.expect.con(model.ans, x, regr.par, 
                      CES = CES, increase = increase, fct1 = fct1, 
                      fct2 = fct2, cc.inf = cc.inf)
                    count <- count + 1
                  }
                }
            }, {
                cc <- ytmp[top]/aa
                if (cont) {
                  if (CES < 0 & cc > 1 + CES) {
                    cc <- 1 + CES - (1 + CES)/100
                  }
                  if (CES > 0 & cc < 1 + CES) {
                    cc <- 1 + CES + CES/100
                  }
                }
                bb <- (aa * (cc - 1) * mean(xtmp))/(mean(ytmp) - 
                  aa) - mean(xtmp)
                if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                if (cc == 0) cc <- 0.1
                if (cc == 1) cc <- 0.9
                CED <- f.inv.con(model.ans = 19, c(aa, bb, cc), 
                  CES)
                if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                  min(x))/3
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc))
            }, {
                dd <- 1
                if (incr.decr.no > 0) CED <- mean(xtmp) else {
                  cc <- ytmp[top]/aa
                  dd <- 1
                  if (cont) {
                    if (CES < 0 & cc > 1 + CES) {
                      cc <- 1 + CES - (1 + CES)/100
                    }
                    if (CES > 0 & cc < 1 + CES) {
                      cc <- 1 + CES + CES/100
                    }
                  }
                  if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                  if (cc == 0) cc <- 0.1
                  if (cc == 1) cc <- 0.9
                  CED <- f.inv.con(model.ans = 20, c(aa, bb, 
                    cc, dd), CES)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                    min(x))/3
                }
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                dt.fr <- data.frame(yyy = ytmp, dose.tmp = xtmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                cc <- res.lm[[1]][1]
                aa <- res.lm[[1]][2]
                bb <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc))
                if (increase == 1) {
                  ans.all$lower <- c(lower.var, eps, eps, eps)
                  ans.all$upper <- c(upper.var, Inf, Inf, Inf)
                }
                if (increase == -1) {
                  ans.all$lower <- c(lower.var, -Inf, eps, eps)
                  ans.all$upper <- c(upper.var, eps, Inf, Inf)
                }
            }, {
                aa <- 1
                bb <- -(logb(aa) - logb(ytmp[top]))/(xtmp[top] - 
                  xtmp[1])
                bb <- 1/bb
                cc <- ytmp[1]
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  cc)
                ans.all$lower <- c(lower.var, eps, -Inf, eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf)
            }, {
                cc <- ytmp[1]
                aa <- xtmp[2]
                bb <- (ytmp[3] - cc)/(xtmp[3] - aa)
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), dd)
                ans.all$lower <- c(lower.var, eps, -Inf, eps, 
                  eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  100)
            }, {
                aa <- ytmp[1]
                cc <- ytmp[top]
                dd <- 1
                bb <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), dd)
                ans.all$lower <- c(lower.var, eps, eps, eps, 
                  eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  100)
            }, {
                dt.fr <- data.frame(yyy = ytmp, dose.tmp = xtmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- res.lm[[1]][1]
                bb <- res.lm[[1]][2]
                cc <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc))
                if (increase == 1) {
                  ans.all$lower <- c(lower.var, eps, eps, eps)
                  ans.all$upper <- c(upper.var, Inf, Inf, Inf)
                }
                if (increase == -1) {
                  ans.all$lower <- c(lower.var, eps, -Inf, eps)
                  ans.all$upper <- c(upper.var, Inf, eps, Inf)
                }
            }, {
                nr.doses <- length(ytmp)
                aa <- ytmp[1]
                bb1 <- mean(xtmp[1:trunc(nr.doses/2)])
                bb2 <- mean(xtmp[trunc(nr.doses/2):nr.doses])
                dd <- 1
                cc1 <- ytmp[nr.doses]/ytmp[1]
                cc2 <- ytmp[nr.doses]/ytmp[1]
                regr.par <- c(rep(aa, nr.aa), bb1, rep(bb2, nr.bb), 
                  rep(2, nr.cc), 0.5, dd)
                ans.all$lower <- c(lower.var, eps, eps, eps, 
                  1, eps, dd.low)
                ans.all$upper <- c(Inf, Inf, Inf, Inf, Inf, 1, 
                  dd.upp)
            }, {
                nr.doses <- length(ytmp)
                aa <- ytmp[1]
                bb1 <- mean(xtmp[1:trunc(nr.doses/2)])
                bb2 <- mean(xtmp[trunc(nr.doses/2):nr.doses])
                dd1 <- 1
                dd2 <- 1
                cc <- ytmp[nr.doses]/ytmp[1]
                regr.par <- c(rep(aa, nr.aa), rep(bb1, nr.bb), 
                  bb2, rep(cc, nr.cc), dd1, dd2)
                ans.all$lower <- c(lower.var, eps, eps, eps, 
                  eps, eps, eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  Inf, 100, 100)
            }, {
                nr.doses <- length(ytmp)
                aa <- ytmp[1]
                bb1 <- mean(xtmp[1:trunc(nr.doses/2)])
                bb2 <- mean(xtmp[trunc(nr.doses/2):nr.doses])
                dd1 <- 1
                dd2 <- 1
                regr.par <- c(rep(aa, nr.aa), bb1, rep(bb2, nr.bb), 
                  rep(2, nr.cc), 0.5, dd1, dd2)
                ans.all$lower <- c(lower.var, eps, eps, eps, 
                  1, eps, dd.low, dd.low)
                ans.all$upper <- c(Inf, Inf, Inf, Inf, Inf, 1, 
                  dd.upp, dd.upp)
            }, {
                aa <- ytmp[1]
                bb <- (ytmp[3] - aa)/(xtmp[3])
                cc <- 0
                dd <- 0
                if (!increase) bb <- -bb
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), dd)
                ans.all$lower <- c(lower.var, -Inf, -Inf, -Inf, 
                  -Inf)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  Inf)
            }, {
                aa <- ytmp[1] + 0.001
                bb <- 1.001
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
                ans.all$lower <- c(lower.var, eps, -Inf)
                ans.all$upper <- c(upper.var, Inf, Inf)
            }, {
                aa <- ytmp[1] + 0.001
                bb <- 1.001
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
                ans.all$lower <- c(lower.var, eps, -Inf)
                ans.all$upper <- c(upper.var, Inf, Inf)
            }, {
                bb <- ytmp[1] + 0.001
                dd <- 1
                aa <- 0
                cc <- 1
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), dd)
                ans.all$lower <- c(lower.var, eps, -Inf, eps, 
                  eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  Inf)
            }, cat(""), cat(""), cat(""), {
                bb <- 1
                aa <- (ytmp[5] - aa)/(xtmp[5] - xtmp[1])
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
                ans.all$lower <- c(lower.var, eps, eps)
                ans.all$upper <- c(upper.var, Inf, Inf)
            }, {
                nr.doses <- length(ytmp)
                aa <- ytmp[1]
                b1 <- mean(xtmp[1:trunc(nr.doses/2)])
                b2 <- b1
                dd <- 1
                cc <- 4
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  b2, rep(cc, nr.cc), dd)
                ans.all$lower <- c(lower.var, eps, eps, eps, 
                  1, eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  Inf, Inf)
            }, {
                dd <- 1
                cc <- ytmp[5] - ytmp[1]
                if (increase == 1) cc <- abs(cc)
                if (increase == -1) cc <- -abs(cc)
                bb <- 1000 * bb
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
                if (increase == 1) {
                  ans.all$lower <- c(lower.var, eps, eps, 0, 
                    eps)
                  ans.all$upper <- c(upper.var, Inf, Inf, 1e+12, 
                    100)
                }
                if (increase == -1) {
                  ans.all$lower <- c(lower.var, eps, eps, -1e+12, 
                    eps)
                  ans.all$upper <- c(upper.var, Inf, Inf, 0, 
                    100)
                }
            }, {
                if (fitted == T) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                }
                if (fitted == F || is.na(bb) || bb == 0) {
                  aa <- mean(logb(yy))
                  bb <- logb(ytmp[5]/ytmp[1])/(xtmp[top] - xtmp[1])
                }
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
                ans.all$lower <- c(lower.var, -Inf, -Inf)
                ans.all$upper <- c(upper.var, Inf, Inf)
            }, {
                nr.doses <- length(ytmp)
                aa <- ytmp[1]
                bb <- mean(xtmp[1:trunc(nr.doses/2)])
                dd <- 1
                cc1 <- ytmp[nr.doses]/ytmp[1]
                cc2 <- ytmp[nr.doses]/ytmp[1]
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  cc1, cc2, dd)
                ans.all$lower <- c(lower.var, eps, eps, eps, 
                  eps, dd.low)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  Inf, dd.upp)
            }, {
                dd <- 1
                if (incr.decr.no > 0) CED <- mean(xtmp) else {
                  cc <- ytmp[top]/aa
                  if (cont) {
                    if ((CES < 0) && (cc > 1 + CES)) {
                      cc <- 1 + CES - (1 + CES)/100
                    }
                    if ((CES > 0) & (cc < 1 + CES)) {
                      cc <- 1 + CES + CES/100
                    }
                  }
                  if (cc == 0) cc <- 0.1
                  if (cc == 1) cc <- 0.9
                  if (!cont) cc <- 0.1
                  CED <- (-(1/bb) * logb((-CES + 1 - cc)/(1 - 
                    cc)))^(1/dd)
                  CED <- abs(CED)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                    min(x))/3
                  if (is.na(CED) || CED == 0) CED <- mean(xtmp)
                }
                CED.vec <- rep(1, nr.bb)
                CED.vec[ref.lev] <- CED
                regr.par <- c(rep(aa, nr.aa), CED.vec, rep(cc, 
                  nr.cc), rep(dd, nr.dd))
            }, {
                dd <- 1
                qq <- 6
                GCED <- mean(x)
                var <- 0.05
                par.start <- c(rep(var, nr.var), rep(aa, nr.aa), 
                  GCED, qq, dd)
                regr.par <- par.start
            }, {
                dd <- 1
                cc <- ytmp[top]/aa
                bb <- xtmp[5]/sf.x
                aa.t <- bb
                bb.t <- 1
                dd.t <- 1
                regr.par <- c(rep(aa, nr.aa), rep(cc, nr.cc), 
                  rep(dd, nr.dd), aa.t, bb.t, dd.t)
            }, {
                bb <- 1/bb
                dd <- 1
                bb <- bb/(log(2))^(1/dd)
                if (fitted == F) {
                  cc <- 1.01 * ytmp[top]/aa
                  if (!is.finite(bb)) bb <- xtmp[5]/100
                } else {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  cc <- ytmp[top]/aa
                }
                bb <- median(x)
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                dd <- 1
                cc <- ytmp[top]/aa
                bb <- xtmp[5]/sf.x
                aa.t <- bb
                bb.t <- 1
                dd.t <- 1
                regr.par <- c(rep(aa, nr.aa), rep(cc, nr.cc), 
                  rep(dd, nr.dd), aa.t, bb.t, dd.t)
            }, {
                dd <- 0.1
                CED <- mean(xtmp)
                if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                  min(x))/3
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(dd, nr.dd))
            }, {
                dd <- 0.1
                if (incr.decr.no > 0) CED <- mean(xtmp) else {
                  cc <- ytmp[top]/aa
                  if (cont && incr.decr.no == 0) {
                    if ((CES < 0) && (cc > 1 + CES)) {
                      cc <- 1 + CES - (1 + CES)/100
                    }
                    if ((CES > 0) & (cc < 1 + CES)) {
                      cc <- 1 + CES + CES/100
                    }
                  }
                  if (cc == 0) cc <- 0.1
                  if (cc == 1) cc <- 0.9
                  if (!cont) cc <- 0.1
                  CED <- mean(xtmp)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                    min(x))/3
                }
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                dd <- 0.1
                CED <- mean(xtmp)
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(dd, nr.dd))
            }, {
                dd <- 0.1
                if (incr.decr.no > 0) CED <- mean(xtmp) else {
                  cc <- ytmp[top]/aa
                  if (cont && incr.decr.no == 0) {
                    if ((CES < 0) && (cc > 1 + CES)) {
                      cc <- 1 + CES - (1 + CES)/100
                    }
                    if ((CES > 0) & (cc < 1 + CES)) {
                      cc <- 1 + CES + CES/100
                    }
                  }
                  if (increase == 1 && cc < 1) cc <- 2
                  if (increase == -1 && cc > 1) cc <- 0.5
                  if (!cont) cc <- 0.1
                  CED <- mean(xtmp)
                }
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
            }, {
                aa <- 1
                bb <- 1
                cc <- 1000
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc))
                ans.all$lower <- c(lower.var, 0, 0, 0)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf)
            }, {
                aa <- 1
                bb <- 1
                cc <- 1000
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc))
                ans.all$lower <- c(lower.var, 0, 0, 0)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf)
            }, {
                if (min(x) != 0) {
                  cat("\n\nATTENTION: lowest value of x (milking) should be zero\n\n")
                  f.press.key.to.continue()
                  return(ans.all)
                }
                if (nr.aa > 1) dum.aa <- yy[x == 0] else dum.aa <- ytmp[1]
                nr.doses <- length(ytmp)
                bb1 <- mean(xtmp[1:trunc(nr.doses/2)])
                bb2 <- mean(xtmp[trunc(nr.doses/2):nr.doses])
                dd1 <- 1
                dd2 <- 1
                nr.b1 <- nr.dd
                nr.b2 <- nr.bb
                regr.par <- c(dum.aa, rep(bb1, nr.b1), rep(bb2, 
                  nr.b2), rep(2, nr.cc), dd1, dd2)
                ans.all$lower <- c(lower.var, eps, 0.001, 0.001, 
                  1.5, dd.low, dd.low)
                ans.all$upper <- c(Inf, Inf, 1000, 1000, 10000, 
                  dd.upp, dd.upp)
            }, {
                FD.V <- 10
                ka <- 0.04
                ke <- 0.01
                regr.par <- c(FD.V, ka, ke)
                ans.all$lower <- c(lower.var, eps, eps, eps)
                ans.all$upper <- c(Inf, Inf, 100, 100)
            }, {
                C1 <- 1
                C2 <- 2
                ka <- 0.04
                lab1 <- 0.04
                lab2 <- 0.01
                regr.par <- c(C1, C2, ka, lab1, lab2)
                ans.all$lower <- c(lower.var, eps, eps, eps, 
                  eps, eps)
                ans.all$upper <- c(Inf, Inf, Inf, 100, 100, 100)
            }, {
                cc <- ytmp[top]/aa
                if (!is.finite(bb)) bb <- xtmp[5]/100
                dd <- 1
                CED <- mean(x)
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
                ans.all$lower <- c(lower.var, -Inf, -Inf, -Inf, 
                  eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  100)
            }, {
                cc <- ytmp[top]/aa
                CED <- mean(x)
                dd <- 1
                aa <- log(aa)
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
                ans.all$lower <- c(lower.var, -Inf, -Inf, -Inf, 
                  eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  100)
            }, {
                cc <- ytmp[top]/aa
                if (!is.finite(bb)) bb <- xtmp[5]/100
                dd <- 1
                aa.can <- aa
                cc.can <- cc
                aa <- log(aa.can * cc.can)
                cc <- qnorm(log(aa.can)/log(aa.can * cc.can))
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
                ans.all$lower <- c(lower.var, -Inf, -Inf, -Inf, 
                  eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  100)
            }, {
                cc <- ytmp[top]/aa
                if (!is.finite(bb)) bb <- xtmp[5]/100
                dd <- 1
                CED <- mean(x)
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
                ans.all$lower <- c(lower.var, -Inf, -Inf, -Inf, 
                  eps)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf, 
                  100)
            }, {
                D <- 1
                V <- 1/max(ytmp)
                lab.min <- -0.5
                lab.plus <- -0.1
                qk <- -2
                regr.par <- c(D, rep(V, nr.aa), rep(lab.min, 
                  nr.bb), rep(lab.plus, nr.cc), rep(qk, nr.dd))
                ans.all$lower <- c(lower.var, 0, 0, -Inf, -Inf, 
                  -Inf)
                ans.all$upper <- c(upper.var, Inf, Inf, 0, 0, 
                  0)
            }, {
                cc <- 10 * max(ytmp)
                k1 <- 0.1
                ke <- 0.2
                regr.par <- c(rep(cc, nr.aa), rep(k1, nr.bb), 
                  rep(ke, nr.dd))
                ans.all$lower <- c(lower.var, 0, 0, 0)
                ans.all$upper <- c(upper.var, Inf, Inf, Inf)
            }, {
                dd <- 1
                if (fitted == F) {
                  cc <- ytmp[top]/aa
                  if (!is.finite(bb)) bb <- xtmp[5]/100
                } else {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  cc <- ytmp[top]/aa
                }
                regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                  rep(cc, nr.cc), rep(dd, nr.dd))
                if (cc < 1) increase <- -1 else increase <- 1
                incr.tmp <- increase + 2
                switch(incr.tmp, {
                  dum.1 <- eps
                  dum.2 <- 1
                }, {
                  dum.1 <- eps
                  dum.2 <- cc.upp
                }, {
                  dum.1 <- 1
                  dum.2 <- cc.upp
                })
                ans.all$lower <- c(lower.var, eps, eps, dum.1, 
                  dd.low)
                ans.all$upper <- c(upper.var, Inf, Inf, dum.2, 
                  dd.upp)
            })
            if (!(model.ans == 11 && model.type == 1)) {
                if (length(xans) > 1) {
                  RPF.vec <- rep(1, (nr.dosecol - 1))
                  regr.par <- c(regr.par, RPF.vec)
                  if (model.ans != 42) 
                    x <- Mx %*% c(1, RPF.vec)
                  ans.all$x <- x
                  ans.all$xy.lim[2:3] <- c(min(ans.all$x), max(ans.all$x))
                }
            }
            ans.all$regr.par <- regr.par
            if (!cont) {
                if (exists("track")) 
                  print("f.start.con:  END.sub")
                ans.all$regr.par <- regr.par
                ans.all$top <- top
                return(ans.all)
            }
            text.par <- f.text.par(ans.all)
            ans.all$nrp <- length(text.par) - nr.var
            if (!(dtype %in% c(10, 15, 250, 260))) {
                expect.trans <- f.expect.con(model.ans, x, regr.par, 
                  fct1, fct2, fct3, fct4, fct5, name = F, CES = CES, 
                  twice = twice, ttt = ttt, y = yy, yy == yy, 
                  increase = increase, x.mn = x.mn, ref.lev = ref.lev, 
                  trace.expect = F, par.start = par.start, sign.q = sign.q, 
                  ans.m6.sd = ans.m6.sd, cc.inf = cc.inf, incr.decr.no = incr.decr.no)
                yy.trans <- yy
                if (sum(expect.trans <= 0) > 0 && !(dtype %in% 
                  c(25, 250))) {
                  cat("\nATTENTION:  model predicts non-negative values \n")
                  f.press.key.to.continue()
                }
                if (dtype == 26) {
                  expect.trans <- sqrt(expect.trans)
                  yy.trans <- sqrt(yy)
                }
                if (dtype %in% c(1, 5, 15)) {
                  expect.trans <- logb(expect.trans)
                  yy.trans <- logb(yy)
                }
                if (sum(is.na(expect.trans)) > 0) {
                  ans.all$par.start <- c(rep(NA, nr.var), regr.par)
                  if (quick.ans > 1) {
                    cat("\n there may be a fitting problem, try single model\n")
                    ans.all$par.start[1] <- 0
                    f.press.key.to.continue()
                  }
                  else if (!WAPP) {
                    cat("\n adjust start values\n")
                    ans.all$text.par <- text.par
                    ans.all$adjust.start <- T
                  }
                  return(ans.all)
                }
                var.start <- var(yy.trans - expect.trans)
            }
            else var.start <- mean(sd2.log)
            if (is.na(var.start)) 
                var.start <- eval(parse(prompt = paste("\ngive value for start value of var (e.g. 0.5)\n")))
            if (model.ans != 47) 
                par.start <- c(rep(var.start, nr.var), regr.par)
            if (ans.all$model.ans == 16) {
                par.start <- MLE
                par.start[length(MLE)] <- BMDratio
                ans.all$regr.par <- par.start[-(1:nr.var)]
                ans.all$nrp <- length(MLE) - nr.var
            }
            ans.all$par.start <- par.start
            ans.all$text.par <- text.par
            ans.all$text.par.0 <- text.par
            ans.all$increase <- increase
            if (plot.type > 0) 
                ans.all$xy.lim <- xy.lim
            if (model.ans == 42) 
                ans.all$plot.type <- 4
            if (plot.type > 8) 
                ans.all$plot.type <- 1
            if (!tmp.quick && plot.type > 0) {
                ans.all$l.ty <- 2
                if (!display_plots == TRUE) {
                  f.graph.window(1)
                }
                f.plot.all(ans.all, new.window = FALSE, no.show = TRUE, display_plots = display_plots)
            }
            switch(fit.ans, loglik.start <- -f.lik.con(ans.all$par.start, 
                x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
                sd2.log, nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                ttt = ttt, trace = trace.tmp, fct4 = fct4, fct5 = fct5, 
                cens.up = cens.up, par.tmp = NA, increase = increase, 
                x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                ans.m6.sd = ans.m6.sd, Mx = Mx, x1 = 1, x2 = x2, 
                cc.inf = cc.inf, incr.decr.no = incr.decr.no), 
                loglik.start <- f.lik.biv(par.start, x, yy, fct1, 
                  fct2, model.ans, trace = F, twice = twice))
            if (is.na(loglik.start) | loglik.start < -1e+20) {
                cat("\nATTENTION: You need to adjust start values before fitting the model...\n")
                ans.all$adjust.start <- T
                if (!WAPP) {
                  f.press.key.to.continue()
                  return(ans.all)
                }
            }
            if (!tmp.quick) {
                if (model.ans == 11) {
                  cat("\n\napproximate value of log-likelihood for full model: ", 
                    round(loglik.start, 2), "\n")
                  cat("\nfit the model for a more precise value of the log-likelihood\n\n")
                }
                else cat("\n\nlog-likelihood at start values: ", 
                  round(loglik.start, 2), "\n\n")
            }
            ans.all$loglik.start <- loglik.start
            ans.all <- f.constr.con(ans.all, tmp.quick = T)
        }
        if (adjust == T) {
            ans.all.plt <- ans.all
            switch(fit.ans, ans.all.plt$regr.par <- par.start[-(1:ans.all$nr.var)], 
                ans.all.plt$regr.par <- par.start)
            nrp <- length(ans.all.plt$regr.par)
            ans <- 1
            if (fitted) 
                loglik.old <- loglik
            else loglik.old <- loglik.start
            ans.all$fitted <- FALSE
            if (sum(is.na(par.start)) == 0) {
                if (ans.plt == 1) {
                  ans.all.plt$heading <- ""
                  ans.all.plt$ans.plt <- 1
                  f.plot.sep(ans.all.plt)
                }
            }
            while (ans[1] != (nrp + (nr.var + 1))) {
                cat("\n=======================================================\n")
                cat(paste("model: ", model.names[model.ans]))
                cat("\n\nold log-likelihood: ", round(loglik.old, 
                  3))
                loglik.new <- NA
                if (sum(is.finite(par.start)) != length(par.start)) 
                  cat("\n\none of the start values has no finite value, first fix that\n\n")
                else switch(fit.ans, loglik.new <- -f.lik.con(par.start, 
                  x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
                  sd2.log, nn, Vdetlim = Vdetlim, CES = CES, 
                  twice = twice, ttt = ttt, trace = trace.tmp, 
                  fct4 = fct4, fct5 = fct5, cens.up = cens.up, 
                  par.tmp = NA, increase = increase, x.mn = x.mn, 
                  ref.lev = ref.lev, sign.q = sign.q, Mx = Mx, 
                  x1 = x1, x2 = x2, ans.m6.sd = ans.m6.sd, cc.inf = cc.inf, 
                  incr.decr.no = incr.decr.no), loglik.new <- f.lik.biv(par.start, 
                  x, yy, fct1, fct2, model.ans, trace = F, twice = twice))
                cat("\nnew log-likelihood: ", round(loglik.new, 
                  3), "\n")
                cat("=======================================================\n\n")
                space <- rep("    current value:      ", (nrp + 
                  nr.var))
                if (model.ans == 11) 
                  cat("\nATTENTION: log-likelihood not implemented for full model\n")
                par.list <- "start values"
                for (ii in 1:length(text.par)) {
                  par.list <- paste(par.list, "\n", ii, "  ", 
                    text.par[ii], " ", signif(par.start[ii], 
                      3))
                }
                par.list <- paste(par.list, "\n", length(text.par) + 
                  1, "  ", "Continue\n")
                cat(par.list)
                ans <- eval(parse(prompt = paste("Which parameter(s) do you want to change?  > ")))
                if (ans[1] != (nrp + nr.var + 1)) {
                  start.tmp <- eval(parse(prompt = paste("Please, give new value for ", 
                    text.par[ans], "  >")))
                  par.start[ans] <- start.tmp
                  if (model.ans == 46) {
                    CED.vec <- par.start[(nr.var + nr.aa + 1):(nr.var + 
                      nr.aa + nr.bb)]
                    CED.ref <- CED.vec[ref.lev]
                    ans.all$CED.ref <- CED.ref
                    ans.all.plt$CED.ref <- CED.ref
                  }
                }
                if (incr.decr.no == 0) 
                  if (cont && model.ans %in% c(6, 14, 15, 19, 
                    20, 24, 25, 52, 54)) {
                    cc <- par.start[nr.var + nr.aa + nr.bb + 
                      1]
                    if (CES < 0 && cc > 1 + CES) {
                      cc <- 1 + CES - (1 + CES)/100
                      cat("\nATTENTION: the value of c you have chosen contradicts value of CES/n")
                    }
                    if (CES > 0 && cc < 1 + CES) {
                      cc <- 1 + CES + CES/100
                      cat("\nATTENTION: the value of c you have chosen contradicts value of CES/n")
                    }
                    par.start[nr.var + nr.aa + nr.bb + 1] <- cc
                  }
                ans.all.plt$par.start <- par.start
                if (length(xans) > 1) {
                  RPF.vec <- c(par.start[(length(par.start) - 
                    nr.dosecol + 2):length(par.start)])
                  ans.all.plt$x <- Mx %*% c(1, RPF.vec)
                  ans.all.plt$xy.lim[2:3] <- c(min(ans.all$x), 
                    max(ans.all$x))
                }
                switch(fit.ans, ans.all.plt$regr.par <- par.start[-(1:nr.var)], 
                  ans.all.plt$regr.par <- par.start)
                if (ans.all.plt$plot.type == 0) 
                  ans.all.plt$plot.type <- 1
                ans.all.plt$show <- ""
                ans.all.plt$l.ty <- 2
                if (ans.plt == 1) 
                  f.plot.sep(ans.all.plt)
                else f.plot.all(ans.all.plt, new.window = FALSE, 
                  no.show = TRUE)
                loglik.old <- loglik.new
            }
            ans.all$par.start <- par.start
            ans.all$adjust.start <- F
            ans.all$loglik.start <- loglik.new
        }
        ans.all$nr.var <- nr.var
        ans.all$npar <- length(ans.all$par.start)
        ans.all$CED <- NA
        ans.all$l.ty <- 1
        if (exists("track")) 
            print("f.start.con:  END")
        return(ans.all)
    })
}



f.start.lm.con <- function(xtmp, ytmp, model.ans, dtype) {
    if (exists("track")) 
        print("f.start.lm.con")
    if (dtype %in% c(25, 250)) 
        dt.fr <- data.frame(yyy = ytmp, dose.tmp = xtmp)
    else if (dtype %in% c(26, 260)) 
        dt.fr <- data.frame(yyy = sqrt(ytmp), dose.tmp = xtmp)
    else dt.fr <- data.frame(yyy = log(ytmp), dose.tmp = xtmp)
    if (!model.ans %in% c(17:25)) {
        res.lm <- lm(yyy ~ dose.tmp, dt.fr)
        aa <- as.numeric(res.lm[[1]][1])
        bb <- as.numeric(res.lm[[1]][2])
        if (dtype %in% c(1, 5, 10, 15)) 
            aa <- exp(aa)
        if (dtype %in% c(26, 260)) 
            aa <- aa^2
    }
    if (model.ans %in% c(17:25)) {
        dt.fr <- data.frame(yyy = 1/ytmp, dose.tmp = xtmp)
        res.lm <- lm(yyy ~ dose.tmp, dt.fr)
        aa <- as.numeric(1/res.lm[[1]][1])
        bb <- as.numeric(1/res.lm[[1]][2]/aa)
    }
    if (bb < 0) 
        increase <- -1
    else increase <- 1
    if (model.ans %in% c(17:20, 22:25)) {
        if (bb < 0) 
            increase <- 1
        else increase <- -1
    }
    if (exists("track")) 
        print("f.start.lm.con:  END")
    return(c(aa, bb, increase))
}



f.report.pars <- function(ans.all) {
    if (exists("track")) 
        print("f.report.pars")
    with(ans.all, {
        if (dtype == 4 && length(unique(fct3)) > 1) 
            nth <- 0
        report.pars <- paste(" estimate for ", text.par, ": ", 
            signif(MLE, 4), "\n")
        if (quick.ans == 1) {
            report.pars <- c(report.pars, "\nlog-likelihood: ", 
                round(ans.all$loglik, 3))
            report.pars <- c(report.pars, "\nAIC: ", round(ans.all$aic, 
                3))
            if (cont && CES != 0) 
                report.pars <- c(report.pars, "\nCES: ", ans.all$CES[1])
            if (!cont && BMR != 0) 
                report.pars <- c(report.pars, "\nBMR: ", ans.all$BMR)
            if (cont && fit.ans == 1) {
                if (max(fct3) == 1) 
                  report.pars <- c(report.pars, "\nestimate for residual C.V.", 
                    fct3.txt[1], ":", paste(round(100 * sqrt(exp(MLE[1]) - 
                      1), 2), "%", sep = ""))
                else for (jj in 1:max(fct3)) report.pars <- c(report.pars, 
                  "\nestimate for residual C.V. in group", fct3.txt[jj], 
                  ":", round(sqrt(exp(MLE[jj]) - 1), 2))
            }
            if (model.type == 2 && !cont) {
                if (dtype == 6) 
                  nrp <- nrp + max(fct3)
                if (lb[nrp + 1] == ub[nrp + 1]) 
                  report.pars[nrp + 1] <- paste("fixed value for th-1:", 
                    lb[nrp + 1], "\n")
                if (lb[nrp + nth + 1] == ub[nrp + nth + 1]) 
                  report.pars[nrp + nth + 1] <- paste("fixed value for sigma:", 
                    lb[nrp + nth + 1], "\n")
            }
        }
        if (exists("track")) 
            print("f.report.pars:   END")
        return(report.pars)
    })
}



f.plot.con <- function(ans.all, sep = FALSE, display_plots = TRUE, save_plots = FALSE) {
    if (exists("track")) 
        print("f.plot.con")
    if (ans.all$plot.type == 0) {
        cat("\nTECHNICAL NOTE: plot.type = 0 (f.plot.con)\n")
        return(ans.all)
    }
    ans.all <- f.data.plt.con(ans.all, sep = sep)
    if (ans.all$ans.m6.sd == 2) 
        ans.all$gr.txt <- ans.all$gr.txt.saved
    with(ans.all, {
        if (plot.type == 9) {
            qqnorm(regr.resid, plot = FALSE)
            qqnorm(regr.resid, xlab = "theoretical quantiles", 
                ylab = "regression residuals", main = "", cex = cex.1, 
                mkh = 0.03, plot = TRUE)
            xxx <- f.distr(regr.resid, Vdetlim, plot.ans = FALSE, 
                ans.distr = 3)
            mu <- xxx[1]
            sig2 <- xxx[2]
            abline(mu, sqrt(sig2))
            f.dtype.txt(dtype)
            return(ans.all)
        }
        bbb <- (x.lim.plt[2] - x.lim.plt[1])/70
        dfr <- points.plt.lst[[1]]
        if (display_plots == TRUE) {
          plot(dfr$x, dfr$y, xlim = x.lim.plt, ylim = y.lim.plt, 
              xlab = xleg, ylab = yleg, cex = 0.5, col = 1, type = "n")
          title(main = heading)
        }
        nr.gr <- length(gr.txt)
        if (plt.mns < 3) 
            for (ii in 1:nr.gr) {
                dfr <- points.plt.lst[[ii]]
                x.tmp <- dfr$x
                y.tmp <- dfr$y
                if (mark[ii] == 1) 
                  cex.1 <- cex.1 * 2
                if (display_plots == TRUE) {
                  points(x.tmp, y.tmp, col = color[ii], pch = mark[ii], 
                  cex = cex.1)
                }
                if (mark[ii] == 1) 
                  cex.1 <- cex.1/2
            }
        if (plt.mns %in% c(1, 3)) 
            for (ii in 1:nr.gr) {
                dfr <- means.plt.lst[[ii]]
                x.tmp <- dfr$x
                y.tmp <- dfr$y
                type.dum <- "p"
                l.ty <- 1
                if (model.ans == 11) {
                  l.ty <- 2
                  type.dum <- "b"
                }
                if (display_plots == TRUE) {
                  points(x.tmp, y.tmp, col = color[ii], pch = mark[ii], 
                  type = type.dum, cex = cex.2)
                }
                ans.CI <- 2
                if (length(x.tmp) > 1000) {
                  ans.CI <- menu(c("no", "yes"), title = "plotting CIs may take some time (or result in unclear plot), \n                           are you sure?")
                }
                if (CI.plt && ans.CI == 2) {
                  L025 <- dfr$CI.low
                  L975 <- dfr$CI.upp
                  if (!is.na(L025[1])) {
                    for (qq in 1:length(x.tmp)) {
                      lines(rep(x.tmp[qq], 2), c(L025[qq], L975[qq]), 
                        col = color[ii])
                      lines(c(x.tmp[qq] - bbb, x.tmp[qq] + bbb), 
                        rep(L025[qq], 2), col = color[ii])
                      lines(c(x.tmp[qq] - bbb, x.tmp[qq] + bbb), 
                        rep(L975[qq], 2), col = color[ii])
                    }
                  }
                }
            }
        if (plot.type == 10) {
            f.dtype.txt(dtype)
            abline(0, 0, lty = 2)
        }
        if (0) 
            if (plot.type == 11) {
                x.plt <- x
                if (x.leg != "") 
                  xleg <- paste("log10 of normalized ", x.leg, 
                    sep = "")
                else xleg <- ""
                y.plt <- yy
                x.plt <- x.plt[!is.na(y.plt)]
                y.plt <- y.plt[!is.na(y.plt)]
                yleg <- "residual"
                yleg <- "residual"
                x.lim.tmp <- log10(x.lim)
                if (!is.finite(x.lim.tmp)[1]) 
                  x.lim.tmp[1] <- log10(dum.contr)
                if (is.na(xy.lim[4])) 
                  y.lim <- c(min(y.plt, na.rm = T), max(y.plt, 
                    na.rm = T))
                y.lim.tmp <- c(min(y.plt, na.rm = T), max(y.plt, 
                  na.rm = T))
            }
        if (max(Vdetlim, na.rm = T) > 0 & plot.type < 9) {
            abline(detlim.plt.1, 0, lty = 3)
            abline(detlim.plt.2, 0, lty = 3)
        }
        if (!is.na(cens.up) & plot.type < 9) {
            lines(c(min(xy.plt$x.plt, na.rm = T), max(xy.plt$x.plt, 
                na.rm = T)), rep(cens.up, 2), lty = 3)
        }
        if (any(label != 0)) {
            if (plot.type == 10) 
                cat("\nthis option is not possible for residuals\n")
            else {
                cat("\n\nUse LEFT mouse button to click on data points that you want to identify\n\n")
                cat("and RIGHT mouse button to stop\n")
                flush.console()
                identify(xy.plt$x.plt, xy.plt$y.plt, labels = label)
            }
        }
        ans.all$nr.gr <- nr.gr
        if (exists("track")) 
            print("f.plot.con:  END")
        return(ans.all)
    })
}



f.data.plt.con <- function(ans.all, sep = FALSE) {
    if (exists("track")) 
        print("f.data.plt.con")
    ans.all$sep <- sep
    with(ans.all, {
        if (is.na(cex.2)) {
            if (plt.mns == 2) 
                cex.2 <- 1
            else cex.2 <- 1.5
        }
        if (dtype %in% c(10, 15, 250, 260)) 
            plt.mns <- 3
        if (covar.no > 0) 
            nr.bb <- max(covariate)
        fct2.save <- fct2
        fct2.txt.save <- fct2.txt
        if (length(xans) > 1) {
            fct2 <- as.numeric(as.factor(displ.fact))
            nr.bb <- max(fct2)
            fct1 <- rep(1, length(fct2))
            nr.aa <- 1
            fct3 <- rep(1, length(fct2))
            fct4 <- rep(1, length(fct2))
            fct5 <- rep(1, length(fct2))
            if (model.ans == 42) 
                nr.bb <- 2
        }
        else {
            if (covar.no > 0) 
                displ.fact <- covariate
            if (max(as.numeric(as.factor(displ.fact))) > 1) {
                fct2 <- as.numeric(as.factor(displ.fact))
                ans.all$fct2 <- fct2
                ans.all$fct2.txt <- levels(factor(displ.fact))
                nr.bb <- max(fct2)
            }
            if (tans == 0) {
                if (ans.m6.sd != 2) {
                  if (model.ans %in% c(31, 33, 57)) 
                    gr.txt <- f.pars.m31(ans.all)$gr.txt
                  else gr.txt <- f.pars(ans.all)$gr.txt
                  ans.all$gr.txt.saved <- gr.txt
                }
                if (ans.m6.sd == 2) 
                  gr.txt <- ans.all$gr.txt.saved
            }
            if (nr.aa == 1 && nr.bb == 1 && nr.cc == 1 && nr.dd == 
                1 && nr.var > 1) {
                fct2 <- fct3
                gr.txt <- fct3.txt
                nr.bb <- max(fct2)
            }
            if (tans > 0 && sep) 
                gr.txt <- ""
        }
        dum.contr <- xy.lim[1]
        dum.zero.resp <- xy.lim[6]
        x.lim <- xy.lim[2:3]
        y.lim <- xy.lim[4:5]
        if (!is.na(cens.up)) 
            y.lim <- c(xy.lim[4], max(xy.lim[5], 1.5 * cens.up))
        if (1) {
            low.x <- x.lim[1]
            upp.x <- x.lim[2]
            yy <- yy[x >= low.x & x <= upp.x]
            nn <- nn[x >= low.x & x <= upp.x]
            if (length(label) > 1) 
                ans.all$label <- label[x >= low.x & x <= upp.x]
            sd2.log <- sd2.log[x >= low.x & x <= upp.x]
            if (length(fct1) > 1) 
                fct1 <- fct1[x >= low.x & x <= upp.x]
            if (length(fct2) > 1) 
                fct2 <- fct2[x >= low.x & x <= upp.x]
            if (length(fct3) > 1) 
                fct3 <- fct3[x >= low.x & x <= upp.x]
            if (length(fct4) > 1) 
                fct4 <- fct4[x >= low.x & x <= upp.x]
            if (length(fct3) > 1) 
                fct5 <- fct5[x >= low.x & x <= upp.x]
            x <- x[x >= low.x & x <= upp.x]
        }
        mark <- c(1:2, 4:25, 33:500)
        f.means <- function(xx, yy, dtype, plot.type, CI = F, 
            sd2.log, nn) {
            if (exists("track2")) 
                print("f.means within f.data.plt.con")
            mean.y <- yy
            if (dtype %in% c(1, 5)) 
                mean.y <- exp(tapply(logb(yy), xx, mean))
            if (dtype %in% c(26)) 
                mean.y <- (tapply(sqrt(yy), xx, mean))^2
            if (dtype %in% c(25)) 
                mean.y <- tapply(yy, xx, mean)
            if (plot.type %in% c(1, 2, 5, 10, 11, 15)) 
                mean.y.plt <- mean.y
            if (plot.type %in% c(3, 4, 6)) 
                mean.y.plt <- log10(mean.y)
            if (plot.type %in% c(7, 8)) 
                mean.y.plt <- sqrt(mean.y)
            if (!CI) 
                return(mean.y.plt)
            if (CI) {
                f.var <- function(dtype, xx, yy, sd2.log, nn) {
                  if (exists("track2")) 
                    print("f.var within f.data.plt.con")
                  yy <- yy[order(xx)]
                  xx <- sort(xx)
                  if (dtype %in% c(10, 15, 250)) {
                    SS <- sum(sd2.log * (nn - 1))
                    df <- sum(nn - 1)
                    if (is.na(df)) 
                      cat("\nATTENTION: variable nn appears to contain NAs;\nthis might be solved by setting ans.all$covar.no equal to 0\n\n")
                    var.within <- SS/df
                  }
                  else {
                    if (dtype %in% c(1, 5)) 
                      yy <- log10(yy)
                    if (dtype %in% c(26)) 
                      yy <- sqrt(yy)
                    mn <- tapply(yy, xx, mean)
                    nn <- tapply(yy, xx, length)
                    Vmn <- rep(mn[1], nn[1])
                    if (length(mn) > 1) {
                      for (ii in 2:length(mn)) Vmn <- c(Vmn, 
                        rep(mn[ii], nn[ii]))
                    }
                    resid <- yy - Vmn
                    df <- sum(nn - 1)
                    if (df == 0) {
                      cat("\nconf intervals cannot be calculated, df = 0\n")
                      var.within <- NA
                    }
                    else var.within <- (sum(nn) - 1) * var(resid)/df
                  }
                  Vsem <- sqrt(var.within/nn)
                  out.lst <- list(Vsem = Vsem, df = df, var.within = var.within)
                  return(out.lst)
                }
                var.out <- f.var(dtype, xx, yy, sd2.log, nn)
                Vsem <- var.out$Vsem
                df <- var.out$df
                if (df > 0) 
                  Vconf <- qt(0.975, df) * Vsem
                else Vconf <- NA
                if (dtype %in% c(1, 5, 10, 15)) {
                  conf.L <- mean.y/exp(Vconf)
                  conf.U <- mean.y * exp(Vconf)
                }
                if (dtype %in% c(25, 250)) {
                  conf.L <- mean.y - Vconf
                  conf.U <- mean.y + Vconf
                }
                if (dtype %in% c(26)) {
                  conf.L <- (sqrt(mean.y) - Vconf)^2
                  conf.U <- (sqrt(mean.y) + Vconf)^2
                }
                if (plot.type %in% c(3, 4, 6)) {
                  conf.L <- log10(conf.L)
                  conf.U <- log10(conf.U)
                }
                if (plot.type %in% c(7, 8)) {
                  conf.L <- sqrt(conf.L)
                  conf.U <- sqrt(conf.U)
                }
                y.lim.plt <- c(min(conf.L), max(conf.U))
                out.lst <- list(conf.L = conf.L, conf.U = conf.U, 
                  y.lim.plt = y.lim.plt, var.within = var.out$var.within)
                return(out.lst)
            }
        }
        detlim.plt.1 <- NA
        detlim.plt.2 <- NA
        if (max(Vdetlim, na.rm = T) > 0) {
            detlim.plt.1 <- min(Vdetlim[Vdetlim > 0], na.rm = T)
            detlim.plt.2 <- max(Vdetlim[Vdetlim > 0], na.rm = T)
        }
        x.plt <- NA
        y.plt <- NA
        if (plot.type == 0) {
            cat("\nTechnical note: plot.type has value 0, in f.data.plt.con\n")
            return(ans.all)
        }
        switch(plot.type, {
            x.plt <- x
            y.plt <- yy
            xleg <- x.leg
            yleg <- y.leg
            x.lim.plt <- x.lim
            y.lim.plt <- y.lim
        }, {
            x.plt <- x
            x.plt[x.plt == 0] <- dum.contr
            x.plt <- log10(x.plt)
            y.plt <- yy
            if (x.leg != "") xleg <- paste("log10-", x.leg, sep = "") else xleg <- ""
            yleg <- y.leg
            y.lim.plt <- y.lim
            x.lim.plt <- log10(x.lim)
            if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
        }, {
            x.plt <- x
            y.plt <- log10(yy)
            cens.up <- log10(cens.up)
            xleg <- x.leg
            if (y.leg != "") yleg <- paste("log10-", y.leg, sep = "") else yleg <- ""
            if (max(Vdetlim, na.rm = T) > 0) {
                detlim.plt.1 <- log10(detlim.plt.1)
                detlim.plt.2 <- log10(detlim.plt.2)
            }
            x.lim.plt <- x.lim
            y.lim.plt <- log10(y.lim)
            if (!is.finite(y.lim.plt)[1]) y.lim.plt[1] <- log10(dum.zero.resp)
        }, {
            x.plt <- x
            x.plt[x.plt == 0] <- dum.contr
            x.plt <- log10(x.plt)
            y.plt <- log10(yy)
            if (is.na(y.lim[1])) y.lim <- xy.lim[4:5]
            cens.up <- log10(cens.up)
            if (x.leg != "") xleg <- paste("log10-", x.leg, sep = "") else xleg <- ""
            if (y.leg != "") yleg <- paste("log10-", y.leg, sep = "") else yleg <- ""
            if (max(Vdetlim, na.rm = T) > 0) {
                detlim.plt.1 <- log10(detlim.plt.1)
                detlim.plt.2 <- log10(detlim.plt.2)
            }
            y.lim.plt <- log10(y.lim)
            x.lim.plt <- log10(x.lim)
            if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
            if (!is.finite(y.lim.plt)[1]) y.lim.plt[1] <- log10(dum.zero.resp)
        }, {
            y.plt <- yy
            x.plt <- sqrt(x)
            if (x.leg != "") xleg <- paste("sqrt-", x.leg, sep = "") else xleg <- ""
            yleg <- y.leg
            y.lim.plt <- y.lim
            x.lim.plt <- sqrt(x.lim)
        }, {
            x.plt <- sqrt(x)
            y.plt <- log10(yy)
            cens.up <- log10(cens.up)
            if (is.na(y.lim[1])) y.lim <- xy.lim[4:5]
            if (y.leg != "") yleg <- paste("log10-", y.leg, sep = "") else yleg <- ""
            if (x.leg != "") xleg <- paste("sqrt-", x.leg, sep = "") else xleg <- ""
            if (max(Vdetlim, na.rm = T) > 0) {
                detlim.plt.1 <- log10(detlim.plt.1)
                detlim.plt.2 <- log10(detlim.plt.2)
            }
            x.lim.plt <- sqrt(x.lim)
            y.lim.plt <- log10(y.lim)
            if (!is.finite(y.lim.plt)[1]) y.lim.plt[1] <- log10(dum.zero.resp)
        }, {
            x.plt <- x
            y.plt <- sqrt(yy)
            cens.up <- sqrt(cens.up)
            if (is.na(y.lim[1])) y.lim <- xy.lim[4:5]
            xleg <- x.leg
            if (y.leg != "") yleg <- paste("sqrt-", y.leg, sep = "") else yleg <- ""
            if (max(Vdetlim, na.rm = T) > 0) {
                detlim.plt.1 <- sqrt(detlim.plt.1)
                detlim.plt.2 <- sqrt(detlim.plt.2)
            }
            x.lim.plt <- x.lim
            y.lim.plt <- sqrt(y.lim)
        }, {
            x.plt <- x
            x.plt[x.plt == 0] <- dum.contr
            x.plt <- log10(x.plt)
            y.plt <- sqrt(yy)
            cens.up <- sqrt(cens.up)
            if (is.na(y.lim[1])) y.lim <- xy.lim[4:5]
            if (y.leg != "") yleg <- paste("sqrt-", y.leg, sep = "") else yleg <- ""
            if (x.leg != "") xleg <- paste("log10-", x.leg, sep = "") else xleg <- ""
            if (max(Vdetlim, na.rm = T) > 0) {
                detlim.plt.1 <- sqrt(detlim.plt.1)
                detlim.plt.2 <- sqrt(detlim.plt.2)
            }
            y.lim.plt <- sqrt(y.lim)
            x.lim.plt <- log10(x.lim)
            if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
        }, {
            ans.all$y.plt <- regr.resid
            ans.all$y.lim.plt <- c(min(ans.all$y.plt, na.rm = T), 
                max(ans.all$y.plt, na.rm = T))
            ans.all$cex.1 <- cex.1
            if (exists("track")) print("f.data.plt.con:   END")
            return(ans.all)
        }, {
            y.plt <- regr.resid
            x.plt <- x
            x.plt[x.plt == 0] <- dum.contr
            x.plt <- log10(x.plt)
            x.plt <- x.plt[!is.na(y.plt)]
            if (x.leg != "") xleg <- paste("log10-", x.leg, sep = "") else xleg <- ""
            yleg <- "residual"
            x.lim.plt <- log10(x.lim)
            if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
            y.lim.plt <- c(min(y.plt, na.rm = T), max(y.plt, 
                na.rm = T))
            if (plt.mns == 3) plt.mns <- 2
        }, {
            x.plt <- x
            if (x.leg != "") xleg <- paste("log10 of normalized ", 
                x.leg, sep = "") else xleg <- ""
            y.plt <- yy
            x.plt <- x.plt[!is.na(y.plt)]
            y.plt <- y.plt[!is.na(y.plt)]
            yleg <- "residual"
            yleg <- "residual"
            x.lim.plt <- log10(x.lim)
            if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
            if (is.na(xy.lim[4])) y.lim <- c(min(y.plt, na.rm = T), 
                max(y.plt, na.rm = T))
            y.lim.plt <- c(min(y.plt, na.rm = T), max(y.plt, 
                na.rm = T))
        })
        shift.plt <- 0
        shift.tmp <- 0
        if (shift > 0) {
            shift.plt <- (max(x.plt, na.rm = T) - min(x.plt, 
                na.rm = T))/shift
            x.lim.plt[2] <- x.lim.plt[2] + shift.plt
        }
        if (gr.txt[1] == "" || is.na(gr.txt[1])) {
            if (covar.no > 0 && !sep) 
                gr.txt <- fct3.txt
            else if (max(fct1) > 1) 
                gr.txt <- fct1.txt
            else if (max(fct2.save) > 1) 
                gr.txt <- fct2.txt
            else if (max(fct3) > 1) 
                gr.txt <- fct3.txt
        }
        if (CI.plt) 
            for (jj in 1:nr.bb) for (ii in 1:nr.aa) for (kk in 1:nr.cc) for (ll in 1:nr.dd) {
                lst <- fct1 == ii & fct2 == jj & fct4 == kk & 
                  fct5 == ll
                x.part <- x.plt[lst]
                y.part <- y.plt[lst]
                y.tmp <- yy[lst]
                if (plot.type == 11) {
                  x.part <- (x.part/max(x.part, na.rm = T))
                  x.part[x.part == 0] <- dum.contr
                  x.part <- log10(x.part)
                }
                if (length(y.part) > 0) {
                  if (dtype %in% c(10, 15, 250)) {
                    sd2.log.part <- sd2.log[lst]
                    nn.part <- nn[lst]
                  }
                  out.lst <- f.means(x.part, y.tmp, dtype, plot.type, 
                    CI = TRUE, sd2.log = sd2.log.part, nn = nn.part)
                  if (!is.na(out.lst$conf.L[1])) {
                    y.lim.CI <- out.lst$y.lim.plt
                    y.lim.plt[1] <- min(y.lim.CI[1], y.lim.plt[1], 
                      na.rm = T)
                    y.lim.plt[2] <- max(y.lim.CI[2], y.lim.plt[2], 
                      na.rm = T)
                  }
                }
            }
        if (nr.aa > 1 && nr.bb > 1 && sum(fct1 != fct2) == 0) {
            nr.aa <- 1
            fct1 <- rep(1, length(x))
        }
        if (nr.aa > 1 & nr.bb > 1 && sum(fct1 != fct2) == 0) {
            nr.aa <- 1
            fct1 <- rep(1, length(x))
        }
        points.plt.lst <- list()
        means.plt.lst <- list()
        if (gr.txt[1] == "") 
            gr.txt <- "all"
        zz <- 0
        for (jj in 1:nr.bb) for (ii in 1:nr.aa) for (kk in 1:nr.cc) for (ll in 1:nr.dd) {
            lst <- fct1 == ii & fct2 == jj & fct4 == kk & fct5 == 
                ll
            x.part <- x.plt[lst]
            y.part <- y.plt[lst]
            y.tmp <- yy[lst]
            if (plot.type == 11) {
                x.part <- (x.part/max(x.part, na.rm = T))
                x.part[x.part == 0] <- dum.contr
                x.part <- log10(x.part)
            }
            if (length(y.part) > 0) {
                zz <- zz + 1
                if (zz > 500) 
                  cat("\nATTENTION: number of subgroups too large for plotting\n\n")
                x.part <- x.part + shift.tmp
                mean.x <- NA
                if (plt.mns %in% c(1, 3)) {
                  if (dtype %in% c(1, 5, 15, 25, 26)) {
                    mean.x <- as.numeric(tapply(x.part, x.part, 
                      mean))
                    mean.y <- f.means(x.part, y.tmp, dtype, plot.type)
                  }
                  if (dtype %in% c(10, 15, 250)) {
                    mean.x <- x.part
                    mean.y <- y.part
                  }
                  if (CI.plt) {
                    if (dtype %in% c(10, 15, 250)) {
                      sd2.log.part <- sd2.log[lst]
                      nn.part <- nn[lst]
                    }
                    out.lst <- f.means(x.part, y.tmp, dtype, 
                      plot.type, CI = TRUE, sd2.log = sd2.log.part, 
                      nn = nn.part)
                    L025 <- out.lst$conf.L
                    L975 <- out.lst$conf.U
                  }
                }
                if (plt.mns == 2 || !CI.plt) {
                  L025 <- rep(NA, length(mean.x))
                  L975 <- rep(NA, length(mean.x))
                }
                points.plt.lst[[gr.txt[zz]]] <- data.frame(x = x.part, 
                  y = y.part)
                if (plt.mns %in% c(1, 3)) 
                  means.plt.lst[[gr.txt[zz]]] <- data.frame(x = as.numeric(mean.x), 
                    y = as.numeric(mean.y), CI.low = as.numeric(L025), 
                    CI.upp = as.numeric(L975))
                shift.tmp <- shift.tmp + shift.plt
            }
        }
        ans.all$xy.plt <- data.frame(x.plt = x.plt, y.plt = y.plt)
        ans.all$points.plt.lst <- points.plt.lst
        ans.all$means.plt.lst <- means.plt.lst
        ans.all$cex.1 <- cex.1
        ans.all$cex.2 <- cex.2
        ans.all$xleg <- xleg
        ans.all$yleg <- yleg
        ans.all$x.lim.plt <- x.lim.plt
        ans.all$y.lim.plt <- y.lim.plt
        ans.all$detlim.plt.1 <- detlim.plt.1
        ans.all$detlim.plt.2 <- detlim.plt.2
        ans.all$gr.txt <- gr.txt
        ans.all$fct2 <- fct2.save
        ans.all$fct2.txt <- fct2.txt.save
        if (exists("track")) 
            print("f.data.plt.con:  END")
        return(ans.all)
    })
}



f.model.specific.results <- function(ans.all) {
    if (exists("track")) 
        print("f.model.specific.results")
    sublist <- list()
    sublist$model.ans <- ans.all$model.ans
    sublist$model.type <- ans.all$model.type
    sublist$modelname <- ans.all$modelname
    sublist$MLE <- ans.all$MLE
    sublist$report.pars <- ans.all$report.pars
    if (length(ans.all$xans) > 1) {
        RPF.vec <- ans.all$regr.par[(length(ans.all$regr.par) - 
            ans.all$nr.dosecol + 2):length(ans.all$regr.par)]
        sublist$x <- ans.all$Mx %*% c(1, RPF.vec)
        xy.lim <- ans.all$xy.lim
        xy.lim[2:3] <- c(min(sublist$x), max(sublist$x))
        sublist$xy.lim <- xy.lim
    }
    sublist$pooled.sd <- ans.all$pooled.sd
    sublist$par.start <- ans.all$MLE
    sublist$loglik.old <- ans.all$loglik
    sublist$nrp <- ans.all$nrp
    sublist$npar <- ans.all$npar
    sublist$npar.aic <- ans.all$npar.aic
    sublist$text.par <- ans.all$text.par
    sublist$lower <- ans.all$lower
    sublist$upper <- ans.all$upper
    sublist$lb <- ans.all$lb
    sublist$ub <- ans.all$ub
    sublist$loglik <- ans.all$loglik
    sublist$aic <- ans.all$aic
    sublist$converged <- ans.all$converged
    sublist$fct1 <- ans.all$fct1
    sublist$fct2 <- ans.all$fct2
    sublist$fct3 <- ans.all$fct3
    sublist$fct4 <- ans.all$fct4
    sublist$fct5 <- ans.all$fct5
    sublist$nr.aa <- max(ans.all$fct1)
    sublist$nr.bb <- max(ans.all$fct2)
    sublist$nr.var <- max(ans.all$fct3)
    sublist$nr.cc <- max(ans.all$fct4)
    sublist$nr.dd <- max(ans.all$fct5)
    sublist$nr.gr <- ans.all$nr.gr
    sublist$gr.txt <- ans.all$gr.txt
    sublist$factor.name <- ans.all$factor.name
    sublist$show <- ans.all$show
    if (ans.all$dtype %in% c(5, 15)) {
        sublist$inter.var <- ans.all$inter.var
        sublist$intra.var <- ans.all$intra.var
    }
    model.ans <- ans.all$model.ans
    model.type <- ans.all$model.type
    if (!(model.ans == 14 & model.type == 1)) {
        sublist$regr.par <- ans.all$regr.par
        sublist$regr.par.matr <- ans.all$regr.par.matr
        sublist$MLE.all <- ans.all$MLE.all
        if (!ans.all$cont) {
            sublist$th.par <- ans.all$th.par
            sublist$sig.par <- ans.all$sig.par
            sublist$nth <- ans.all$nth
            sublist$sens.lev <- ans.all$sens.lev
            sublist$rank.low <- ans.all$rank.low
        }
        sublist$CED <- ans.all$CED
        sublist$CED.matr <- ans.all$CED.matr
        sublist$CED.boot <- ans.all$CED.boot
        sublist$CES <- ans.all$CES
        sublist$ced.table <- ans.all$ced.table
        sublist$loglik.summary <- ans.all$loglik.summary
        sublist$conf.int <- ans.all$conf.int
        if (ans.all$dtype == 3) 
            sublist$conf.int.matr <- ans.all$conf.int.matr
        if (!ans.all$cont) 
            sublist$pi.bg.matr <- ans.all$pi.bg.matr
        sublist$response.matr <- ans.all$response.matr
        sublist$trend <- ans.all$trend
        sublist$regr.resid <- ans.all$regr.resid
        sublist$lines.plt.lst <- ans.all$lines.plt.lst
        sublist$cedes.plt.lst <- ans.all$cedes.plt.lst
        sublist$report.pars <- ans.all$report.pars
        if (ans.all$dtype %in% c(1, 5)) 
            sublist$outliers <- ans.all$outliers
    }
    if (exists("track")) 
        print("f.model.specific.results:  END")
    return(sublist)
}



f.ced.con <- function(ans.all) {
    if (exists("track")) 
        print("f.ced.con")
    with(ans.all, {
        CED <- numeric()
        if (dtype == 4 && max(fct3) > 1) {
            regr.par.matr <- matrix(regr.par.matr[1, ], nrow = 1)
        }
        if (0) 
            if (length(xans) > 1) 
                regr.par.matr <- f.pars(ans.all)$regr.par.matr
        if (model.ans == 48) {
            ans.all$fct2 <- as.numeric(factor(ttt))
            lst.cxt <- f.model.cxt(ans.all)
            ans.all$regr.par <- lst.cxt$regr.par
            ans.all$model.ans <- lst.cxt$model.ans
            ans.all <- f.pars(ans.all)
            model.ans <- ans.all$model.ans
            regr.par.matr <- ans.all$regr.par.matr
        }
        if (model.ans == 16) 
            CED <- rep(0, nr.bb)
        else if (model.ans == 6 && ans.m6.sd == 2) 
            CED <- f.inv.con(model.ans, regr.par.matr[1, -1], 
                CES, ans.m6.sd = 2)
        else if (model.ans == 47) 
            CED <- regr.par.matr[, 3]
        else if (model.ans == 49) 
            CED <- regr.par.matr[, 2]
        else {
            nr.row <- length(regr.par.matr[, 1])
            for (ii in 1:nr.row) {
                par.tmp <- regr.par.matr[ii, ]
                CED[ii] <- f.inv.con(model.ans, par.tmp, CES, 
                  max.x = max.x, ES.abs = ES.abs, m6.ED50 = m6.ED50)
            }
        }
        CED.lst <- list(CED = CED)
        if (exists("track")) 
            print("f.ced.con: END")
        return(CED.lst)
    })
}



f.inv.con <- function(model.ans, params, CES, max.x = NA, ES.abs = 1, ans.m6.sd = 1, 
    dtype = 1, m6.ED50 = FALSE) {
    if (exists("track")) 
        print("f.inv.con")
    CED <- NA
    CES.tmp <- CES
    aa <- params[1]
    bb <- params[2]
    par3 <- params[3]
    par4 <- params[4]
    if (m6.ED50 && model.ans == 6) {
        CED <- bb * (log(2)^(1/par4))
        return(CED)
    }
    switch(model.ans, {
        CED <- NA
        cat("\nfor null model CED is not defined\n")
    }, {
        if (bb < 0) CES.tmp <- -abs(CES)
        CED <- (1/bb) * logb(CES.tmp + 1)
    }, {
        if (bb < 0) CES.tmp <- -abs(CES)
        CED <- (1/bb) * logb(CES.tmp + 1)
        CED <- CED^(1/par3)
    }, {
        if ((par3 < 1)) CES.tmp <- -abs(CES)
        dum <- (CES.tmp + 1 - par3)/(1 - par3)
        dum[dum < 0] <- NA
        if (dtype == 3) CED <- -(1/bb) * logb(dum) else {
            if (sum(is.na(dum)) > 0) {
                cat("\nATTTENTION: parameter c does not allow chosen value for CES\n")
                CED <- NA
            } else CED <- -(1/bb) * logb(dum)
        }
    }, {
        if ((par3 < 1)) CES.tmp <- -abs(CES)
        dum <- (CES.tmp + 1 - par3)/(1 - par3)
        dum[dum < 0] <- NA
        if (dtype == 3) CED <- (-(1/bb) * logb(dum))^(1/par4) else {
            if (sum(is.na(dum)) > 0) {
                cat("\nATTTENTION: parameter c does not allow chosen value for CES\n")
                CED <- NA
            } else CED <- (-(1/bb) * logb(dum))^(1/par4)
        }
    }, {
        if ((par3 < 1)) CES.tmp <- -abs(CES)
        if (ans.m6.sd == 1) CED <- bb * (-log(1 - log(CES.tmp + 
            1)/log(par3)))^(1/par4)
        if (ans.m6.sd == 2) CED <- bb * (-log(1 - CES))^(1/par4)
        if (0) {
            print("f.inv.con    CES.tmp")
            print(CES.tmp)
            print(CES.tmp + 1)
            print(log(par3))
            print("model 6")
            print(-log(1 - log(CES.tmp + 1)/log(par3)))
            print((1 - log(CES.tmp + 1)/log(par3)))
            print(log(CES.tmp + 1)/log(par3))
            print(CED)
        }
    }, {
        if (bb < 0) CES.tmp <- -abs(CES)
        CED <- bb * (logb(CES.tmp + 1))
    }, {
        if (bb < 0) CES.tmp <- -abs(CES)
        if (bb > 0) CED <- bb * (logb(CES.tmp + 1))^(1/par3) else CED <- -bb * 
            (-logb(CES.tmp + 1))^(1/par3)
    }, {
        if (par3 < 1) CES.tmp <- -abs(CES)
        dum <- (CES.tmp + 1 - par3)/(1 - par3)
        dum[dum < 0] <- NA
        if (dtype == 3) CED <- bb * ((-logb(dum))) else {
            if (sum(is.na(dum)) > 0) {
                cat("\nATTTENTION: parameter c does not allow chosen value for CES\n")
                CED <- NA
            } else CED <- bb * ((-logb(dum)))
        }
    }, {
        if (par3 < 1) CES.tmp <- -abs(CES)
        dum <- (CES.tmp + 1 - par3)/(1 - par3)
        dum[dum < 0] <- NA
        if (dtype == 3) {
            CED <- (-(1/bb) * logb(dum))^(1/par4)
            CED <- bb * ((-logb(dum))^(1/par4))
        } else {
            if (sum(is.na(dum)) > 0) {
                cat("\nATTTENTION: parameter c does not allow chosen value for CES\n")
                CED <- NA
            } else {
                CED <- (-(1/bb) * logb(dum))^(1/par4)
                CED <- bb * ((-logb(dum))^(1/par4))
            }
        }
    }, , CED <- bb, CED <- bb, CED <- bb, CED <- bb, {
    }, {
        if (bb > 0) CES.tmp <- -abs(CES)
        CED <- -bb * (CES.tmp/(CES.tmp + 1))
    }, {
        if (bb > 0) CES.tmp <- -abs(CES)
        CED <- ((-sign(bb) * CES.tmp)/(1 + CES.tmp))^(1/par3) * 
            (abs(bb))
    }, {
        if ((par3 < 1)) CES.tmp <- -abs(CES)
        CED <- sign(bb) * bb * ((CES.tmp)/(par3 - 1 - CES.tmp))
    }, {
        if ((par3 < 1)) CES.tmp <- -abs(CES)
        if (par3 == 0) CED <- ((-(bb^par4) * CES.tmp)/(1 + CES.tmp))^(1/par4) else CED <- bb * 
            ((CES.tmp)/(par3 - 1 - CES.tmp))^(1/par4)
    }, {
        if ((par3 < 1)) CES.tmp <- -abs(CES)
        dum <- log(CES.tmp + 1)/log(par3)
        CED <- bb * (dum/(1 - dum))^(1/par4)
    }, CED <- bb, CED <- bb, CED <- bb, CED <- bb, {
        if (ES.abs == 1) {
            CED <- ((CES.tmp - par3)/aa)^(1/bb)
            CED[CED <= 0] <- 1e-06
        } else {
            if (aa < 0) CES.tmp <- -abs(CES)
            CED <- ((par3 * CES.tmp)/aa)^(1/bb)
        }
    }, cat("calculation of CED not implemented"), cat("calculation of CED not implemented"), 
        cat("calculation of CED not implemented"), {
            if (ES.abs == 1) {
                CED <- ((CES.tmp - aa)/bb)^(1/par3)
                CED[CED <= 0] <- 1e-06
            } else {
                if (bb < 0) CES.tmp <- -abs(CES)
                CED <- ((aa * CES.tmp)/bb)^(1/par3)
            }
        }, {
            f.dum <- function(x, params, CES) {
                b1 <- params[2]
                b2 <- params[3]
                c1 <- params[4]
                c2 <- params[5]
                dd <- params[6]
                y.expect.1 <- c1 - (c1 - 1) * exp(-(x/b1)^dd)
                y.expect.2 <- c2 - (c2 - 1) * exp(-(x/b2)^dd)
                y.expect <- (y.expect.1 * y.expect.2)
                return(y.expect - CES - 1)
            }
            CED <- NA
            uniroot.out <- uniroot(f.dum, interval = c(0, max.x), 
                params = params, CES = CES)
            CED <- uniroot.out$root
        }, cat(""), cat(""), cat(""), CED <- -log(CES + 1)/log(bb), 
        cat(""), cat(""), cat(""), cat(""), cat(""), {
            CED <- (CES.tmp/aa)^(1/bb)
        }, {
            cat("")
        }, {
            dum <- abs(CES.tmp - par3/aa)
            CES.tmp <- abs(CES.tmp)
            CED <- bb * CES.tmp^(1/par4)/dum^(1/par4)
        }, cat(""), cat(""), {
            CED <- bb
        }, {
        }, {
            if ((par3 < 1)) CES.tmp <- -abs(CES)
            CED <- bb * (-log(1 - log(CES.tmp + 1)/log(par3)))^(1/par4)
        })
    if (model.ans %in% 51:54) 
        CED <- bb
    if (model.ans == 60) {
        dum <- log(1 - log(CES + 1)/(aa * (par3 - 1)))
        CED <- (-1/bb * dum)^(1/par4)
    }
    if (model.ans == 61) 
        CED <- bb
    if (model.ans == 63) 
        CED <- bb
    if (model.ans == 62) {
        cc <- par3
        dd <- par4
        dum <- qnorm((1/aa) * log(CES + 1) + pnorm(cc))
        CED <- ((dum - cc)/bb)^(1/dd)
    }
    if (dtype != 3 && !ES.abs) 
        if (is.na(CED) && bb < 1e-10) 
            CED <- 1e+10
    if (exists("track")) 
        print("f.inv.con : END")
    return(CED)
}



f.resid.con <- function(ans.all) {
    if (exists("track")) 
        print("f.resid.con")
    with(ans.all, {
        pred.value <- f.expect.con(model.ans, x, regr.par = regr.par, 
            fct1 = fct1, fct2 = fct2, fct3 = fct3, fct4 = fct4, 
            fct5 = fct5, name = F, CES = CES, twice = twice, 
            ttt = 0, y = yy, increase = increase, ref.lev = ref.lev, 
            x.mn = x.mn, ans.m6.sd = ans.m6.sd, sign.q = sign.q, 
            x1 = x1, x2 = x2, cc.inf = cc.inf)
        if (dtype %in% c(1, 5, 10, 15)) 
            regr.resid <- logb(yy) - logb(pred.value)
        if (dtype %in% c(25, 250)) 
            regr.resid <- pred.value - yy
        if (dtype %in% c(26, 260)) 
            regr.resid <- sqrt(yy) - sqrt(pred.value)
        ans.all$regr.resid.raw <- regr.resid
        var.subgr <- MLE[1:nr.var]
        var.vector <- rep(0, length(x))
        if (nr.var == 1) 
            var.vector <- var.vector + var.subgr[1]
        else for (ii in (1:nr.var)) var.vector <- var.vector + 
            var.subgr[ii] * (fct3 == ii)
        ans.all$regr.resid <- regr.resid/sqrt(var.vector)
        ans.all$pred.value <- pred.value
        if (exists("track")) 
            print("f.resid.con:  END")
        if (quick.ans > 1 && dtype %in% c(1, 5, 25, 26)) 
            ans.all <- f.outliers(ans.all)
        return(ans.all)
    })
}



f.outliers <- function(ans.all) {
    if (exists("track")) 
        print("f.outliers")
    with(ans.all, {
        length.old <- length(regr.resid)
        sig2 <- 1
        mu <- mean(regr.resid)
        limit.tmp <- f.grubb(length.old)
        if (dtype %in% c(10, 15)) {
            cat(" \n\nATTENTION:  identifying outliers is not valid for summary data \n\n")
            f.press.key.to.continue()
        }
        if (output) 
            cat(" \nThe Grubb outlier for this sample size is:", 
                round(limit.tmp, 5), "\n")
        limit.low <- -limit.tmp * sqrt(sig2) + mu
        limit.up <- limit.tmp * sqrt(sig2) + mu
        outlier.vec <- regr.resid[(regr.resid < limit.low) & 
            (y > -1000)]
        outlier.vec <- c(outlier.vec, regr.resid[(regr.resid > 
            limit.up) & (y > -1000)])
        resid.1 <- regr.resid[y > -1000]
        resid.2 <- resid.1[(resid.1 > limit.low) & (resid.1 < 
            limit.up)]
        length.old <- length(resid.1)
        length.new <- length(resid.2)
        if (output) 
            if (length.old == length.new) 
                cat(" \nNo outliers detected \n")
        if (length.old > length.new) {
            lst <- (regr.resid < limit.low) | (regr.resid > limit.up) & 
                (y > -1000)
            if (output && sum(lst) > 0) 
                cat(" \n", sum(lst), "outliers detected, with values for x and y: \n")
            numbering <- 1:length(lst)
            rrr <- numbering[lst]
            for (ii in rrr) cat(paste(x[ii], y[ii]), " \n\n")
            ans.all$out.ans <- 2
            if (quick.ans == 1 && length(outlier.vec) > 0 && 
                length(xans) == 1 && !WAPP) 
                ans.all$out.ans <- menu(c("yes", "no"), title = " \n\nDo you want to remove these outliers and refit the model?")
            if (quick.ans == 1) 
                if (ans.all$out.ans == 2 && ans.all$plot.type %in% 
                  c(9:10)) 
                  title(main = "outliers included")
            if (ans.all$out.ans == 1) {
                cat(" \n\nThe following data are removed as outliers:\n\n")
                cat(paste(x[lst], y[lst]))
                ans.all$y <- y[!lst]
                ans.all$yy <- yy[!lst]
                ans.all$xy.lim[4] <- min(ans.all$yy)
                ans.all$x <- x[!lst]
                ans.all$fct1 <- fct1[!lst]
                ans.all$fct2 <- fct2[!lst]
                ans.all$fct3 <- fct3[!lst]
                fct4 <- fct4[!lst]
                ans.all$fct4 <- fct4[!is.na(fct4)]
                fct5 <- fct5[!lst]
                ans.all$fct5 <- fct5[!is.na(fct5)]
                ans.all$data.0 <- data.0[!lst, ]
                ans.all$regr.resid <- regr.resid[!lst]
                f.press.key.to.continue("\n\nthe model will be refitted  \n\n")
                if (first.outliers) {
                  Moutliers <- data.0[1, ]
                  Moutliers[1, ] <- rep(NA, length(data.0[1, 
                    ]))
                  ans.all$Moutliers <- Moutliers
                  ans.all$first.outliers <- FALSE
                }
                ans.all$Moutliers <- rbind(ans.all$Moutliers, 
                  data.0[lst, ])
            }
        }
        if (exists("track")) 
            print("f.outliers:   END")
        return(ans.all)
    })
}



f.factorname <- function(ans.all) {
    if (exists("track")) 
        print("f.factorname")
    with(ans.all, {
        factor.name <- ""
        fact1.name <- ""
        fact2.name <- ""
        fact3.name <- ""
        fact4.name <- ""
        fact5.name <- ""
        if (quick.ans > 1) {
            fct1.no <- 0
            fct2.no <- 0
            fct3.no <- 0
            if (nr.aa > 1) 
                fct1.no <- covar.no
            if (nr.bb > 1) 
                fct2.no <- covar.no
            if (nr.var > 1) 
                fct3.no <- covar.no
        }
        if (model.ans != 55) {
            if (fct1.no > 0) {
                fact1.name <- paste("\nfact-a: ", varnames[fct1.no])
                factor.name <- paste(factor.name, fact1.name)
            }
            if (fct2.no > 0) {
                fact2.name <- paste("\nfact-b: ", varnames[fct2.no])
                factor.name <- paste(factor.name, fact2.name)
            }
            if (fct3.no > 0) {
                fact3.name <- paste("\nfact-var: ", varnames[fct3.no])
                if (dtype == 4) 
                  fact3.name <- paste("\nfact-theta: ", varnames[fct3.no])
                factor.name <- paste(factor.name, fact3.name)
            }
            if (fct4.no > 0) {
                fact4.name <- paste("\nfact-c: ", varnames[fct4.no])
                factor.name <- paste(factor.name, fact4.name)
            }
            if (fct5.no > 0) {
                fact5.name <- paste("\nfact-d: ", varnames[fct5.no])
                factor.name <- paste(factor.name, fact5.name)
            }
        }
        if (model.ans == 55) {
            if (fct1.no > 0) {
                fact1.name <- paste("\nfact-Km: ", varnames[fct1.no])
                factor.name <- paste(factor.name, fact1.name)
            }
            if (fct2.no > 0) {
                fact2.name <- paste("\nfact-Vmax: ", varnames[fct2.no])
                factor.name <- paste(factor.name, fact2.name)
            }
            if (fct4.no > 0) {
                fact4.name <- paste("\nfact-V: ", varnames[fct4.no])
                factor.name <- paste(factor.name, fact4.name)
            }
        }
        if (exists("track")) 
            print("f.factorname  END")
        return(factor.name)
    })
}



f.hit.constr <- function(ans.all) {
    if (exists("track")) 
        print("f.hit.constr")
    with(ans.all, {
        if (!(model.ans == 11 && model.type == 1) && !(model.ans == 
            14 && model.type == 1)) {
            hit.constr <- FALSE
            if (0) {
                print("f.hit.constr")
                print(model.ans)
                print(cbind(lb, MLE, ub))
            }
            lst <- (MLE == lb) * !(MLE == ub)
            lst <- lst == 1
            if (model.type == 2) 
                lst[length(lst)] <- F
            if (model.type == 2 && dtype %in% c(4, 6)) 
                lst[length(lst) - 1] <- F
            if (any(lst)) {
                hit.constr <- TRUE
                cat("\nthe following parameter estimate was equal to the lower constraint:\n")
                print(c(text.par[lst], signif(MLE[lst], 5)))
            }
            lst <- (MLE == ub) * !(MLE == lb)
            lst <- lst == 1
            if (model.type == 2) 
                lst[length(lst)] <- F
            if (model.type == 2 && dtype %in% c(4, 6)) 
                lst[length(lst) - 1] <- F
            if (any(lst)) {
                hit.constr <- TRUE
                cat("\nthe following parameter estimate was equal to the upper constraint:\n")
                print(c(text.par[lst], signif(MLE[lst], 5)))
            }
        }
        if (exists("track")) 
            print("f.hit.constr:   END")
        return(invisible())
    })
}



f.refit.nes <- function(ans.all, interactive_mode = T, results_env = NULL) {
    if (exists("track")) 
        print("f.refit.nes")
    ans.all <- f.nes(ans.all)
    ans.all <- with(ans.all, {
        if (increase == 1 && CES < 0) {
            cat("\nproblem in sign of CES (f.refit.nes)\n")
            f.press.key.to.continue()
        }
        if (increase == -1 && CES > 0) {
            cat("\nproblem in sign of CES (f.refit.nes)\n")
            f.press.key.to.continue()
        }
        if (increase == 0) 
            cat("\nATTENTION:  variable increase has value 0 (f.refit.nes)\n")
        if (output) 
            cat("\nEndpoint-specific CES: ", round(CES, 3), "\n\n")
        if (1) {
            cc.OK <- f.check.cc(ans.all)
            if (!cc.OK) {
                note.tmp <- paste(y.leg, ": parameter cc too close to CES; no BMD CI calculated")
                ans.all$notes <- paste(ans.all$notes, "\n", note.tmp)
                ans.all$conf.int <- matrix(NA, ncol = 2, nrow = max(fct2))
                return(ans.all)
            }
        }
        ii <- length(Vloglik) + 1
        regr.par <- MLE[-(1:nr.var)]
        ans.all$regr.par <- regr.par
        ans.all$model.type <- 1
        ans.all <- f.pars(ans.all)
        ans.all$max.x <- NA
        CED.lst <- f.ced.con(ans.all)
        text.CED <- CED.lst$text.CED
        CED <- CED.lst$CED
        CED[is.na(CED)] <- mean(x)
        CED[is.infinite(CED)] <- 1e+12
        if (is.na(CED[1])) {
            ans.all$CED <- CED
            cat("\n\nATTENTION:  model could not be refitted, no conf,int derived\n\n")
            return(ans.all)
        }
        CED <- abs(CED)
        ans.all$par.start <- MLE
        ans.all$par.start[(nr.var + nr.aa + 1):(nr.var + nr.aa + 
            nr.bb)] <- CED[1:nr.bb]
        if (0) {
            loglik.try <- NA
            count <- 0
            while (is.na(loglik.try)) {
                count <- count + 1
                loglik.try <- -f.lik.con(ans.all$par.start, x, 
                  y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
                  sd2.log, nn, Vdetlim = Vdetlim, CES = CES, 
                  twice = twice, ttt = ttt, trace = F, fct4 = fct4, 
                  fct5 = fct5, cens.up = cens.up, par.tmp = NA, 
                  increase = increase, ref.lev = ref.lev, sign.q = sign.q, 
                  ans.m6.sd = ans.m6.sd, Mx = Mx, x1 = x1, x2 = x2)
                if (is.na(loglik.try)) {
                  if (any(CED > 100)) 
                    CED[CED > max(x)] <- max(x) * 10
                  else CED <- CED * 10
                  ans.all$par.start[(nr.var + nr.aa + 1):(nr.var + 
                    nr.aa + nr.bb)] <- CED[1:nr.bb]
                }
                if (count > 10) {
                  cat("\n\nATTENTION: confidence interval could not be calculated\n")
                  ans.all$conf.int <- matrix(rep(NA), ncol = 2)
                  return(ans.all)
                }
            }
        }
        if (is.na(CED[1])) 
            no.CI <- T
        ans.all$par.start <- ans.all$par.start * 1.05
        model.txt[ii] <- modelname
        if (model.ans %in% c(15, 25, 52, 54)) {
            cc <- regr.par[nr.aa + nr.bb + 1]
            cc.dum <- abs(log(cc))
            ES.dum <- abs(log(CES + 1))
            if (ES.dum >= cc.dum) {
                cat("\n\n CES exceeds value of c, constraints will be adjusted; note that refitted model may not be optimal\n\n\n")
                ans.all <- f.constr.con(ans.all, tmp.quick = TRUE)
            }
            if (0) 
                if (ES.dum >= cc.dum) {
                  cat("\n\n model cannot be re-fitted, CES exceeds value of c\n\n\n")
                  ans.all$do.MA <- FALSE
                  return(ans.all)
                }
        }
        if (!MA.running) 
            cat("\nre-fitting selected model", model.txt[ii], 
                "with endpoint-specific CES", round(CES, 3), 
                "\n")
        ans.all <- f.qfit(ans.all, plot.type = 0)
        if (interactive_mode == TRUE) {
            assign("last.fit", ans.all, immed = T, pos = 1)
        } else {
            #assign(ans.all$modelname, ans.all, envir = results_env)
        }
        if (ans.all$MLE[1] == 0) {
            cat("\n re-fitting model failed\n")
            cat("CES:", CES, "\n")
            ans.all$fitted <- FALSE
            return(ans.all)
        }
        Vloglik[ii] <- round(ans.all$loglik, 2)
        if (0) 
            if (!gui && length(ans.all$MLE) > 10) {
                f.store.results(ans.all, "refit.tmp")
                cat("\nNOTE: refitted  model was stored in refit.tmp\n")
            }
        if (Vloglik[ii] < loglik - 0.1) {
            cat("\nATTENTION: the log-likelihood for the refitted model has decreased:", 
                Vloglik[ii], ", was:", loglik, "\n")
            cat("\nNew attempt of fitting the model ...\n")
            ans.all$par.start <- NA
            ans.all <- f.qfit(ans.all, plot.type = 0)
            f.hit.constr(ans.all)
            Vnpar[ii] <- ans.all$npar
            Vconverged[ii] <- ans.all$converged
            Vloglik[ii] <- round(ans.all$loglik, 2)
        }
        if (Vloglik[ii] < loglik - 0.1) {
            cat("\nATTENTION: the log-likelihood for the refitted model is still lower,", 
                Vloglik[ii], ", while it was:", loglik, "\n")
            if (model.ans %in% c(14, 15, 24, 25)) {
                cat("\nthis might be due to parameter c being too close to CES\n")
            }
        }
        MLE <- ans.all$MLE
        regr.par <- MLE[-(1:max(fct3))]
        CED <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
        if (nr.bb == 1) 
            CED <- rep(CED, nr.aa)
        if (dtype == 5) {
            if (nest.no != 0) {
                if (max(Vdetlim, na.rm = T) != 0) 
                  cat("\nlitter effects not implemented for nonzero detection limit\n")
                else {
                  ans.all <- f.nested.con(ans.all)
                }
            }
            else dtype <- 1
        }
        ans.all$regr.par <- regr.par
        ans.all$regr.par.matr <- f.pars(ans.all)$regr.par.matr
        Vnpar[ii] <- ans.all$npar
        Vconverged[ii] <- ans.all$converged
        Vloglik[ii] <- round(ans.all$loglik, 2)
        Vaic[ii] <- 2 * Vnpar[ii] - 2 * Vloglik[ii]
        if (!MA.running) 
            cat(paste("\n", model.txt[ii], tb, tb, Vconverged[ii], 
                tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii], 
                "")
        if (interactive_mode == TRUE) {
            assign("last.fit", ans.all, immed = T, pos = 1)
        } else {
            #assign(ans.all$modelname, ans.all, envir = results_env)
        }
        if (!no.CI) {
            if (dtype == 5) {
                ans.all$plot.ans <- 1
                ans.all$CED <- CED
                ans.all <- f.mm7.con(ans.all)
                conf.int <- ans.all$conf.int
            }
            else {
                if (ans.all$loglik == 1e-12) 
                  cat("\nno finite log-likelihood found, CI will not be calculated\n")
                else {
                  ans.all$trace <- F
                  ans.all$trace.plt <- F
                  ans.all <- f.CI(ans.all)
                  conf.int <- ans.all$conf.int
                  if (ans.all$update) 
                    MLE <- ans.all$MLE
                }
            }
            # If you want to do the CES recalculation you need to uncomment these lines
            # Something weird about ans.all at this point though. Probably best to avoid this.
            # ans.all$CED <- CED # This is probably not a good idea, but necessary to get the CED to match the output below, which also matches the figures. No idea where the proper ans.all to capture for this is.
            # assign(ans.all$modelname, ans.all, envir = results_env)
            regr.par <- MLE[-(1:max(fct3))]
            CED <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
            if (nr.bb == 1) 
                CED <- rep(CED, nr.aa)
            if (model.ans > 1) {
                CED.uniq <- unique(CED)
                if (!no.CI && cont) 
                  if (length(CED.uniq) > 1) 
                    for (qq in 1:max(fct2)) cat("\nthe CED (in orig. units) and the", 
                      100 * conf.lev, "% confidence interval for group", 
                      covar.txt[qq], "is: \n", signif(sf.x * 
                        CED[qq], 5), "\n", round(conf.int[qq, 
                        1], 5), "\n", round(conf.int[qq, 2], 
                        5), "\n")
                  else cat("\nthe CED (in orig. units) and the", 
                    100 * conf.lev, "% confidence interval is: \n", 
                    signif(sf.x * CED[1], 5), "\n", round(conf.int[1, 
                      1], 5), "\n", round(conf.int[1, 2], 5), 
                    "\n")
            }   
        }
        ans.all$MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
        ans.all$CED <- CED
        ans.all$regr.par <- regr.par
        ans.all$regr.par.matr <- f.pars(ans.all)$regr.par.matr
        ans.all$ab.txt <- ab.txt
        ans.all$Vloglik <- Vloglik
        ans.all$Vaic <- Vaic
        ans.all$Vnpar <- Vnpar
        ans.all$Vconverged <- Vconverged
        ans.all$MLE <- MLE
        if (!cont & dtype != 3 & model.ans > 1) 
            ans.all$CED <- CED
        if (0) 
            if (!gui && length(ans.all$MLE) > 10) {
                f.store.results(ans.all, "refit.tmp")
                cat("\nNOTE: confidence intervals were stored in refit.tmp\n")
            }
        ans.all$model.txt <- model.txt
        ans.all.tmp <- ans.all
    }, return(ans.all.tmp)
    
)
    if (exists("track")) 
        print("f.refit.nes:  END")
    return(ans.all)
    
}



f.nes <- function(ans.all) {
    if (exists("track")) 
        print("f.nes")
    with(ans.all, {
        var.within.pooled <- f.var.pooled(ans.all)
        sd.0 <- log(1.05)
        pooled.sd <- sqrt(var.within.pooled)
        ES.0 <- 1.05
        ES.new = exp(log(ES.0) * pooled.sd/sd.0)
        if (increase == 1) {
            CES <- ES.new - 1
            ans.all$lb[nr.var + nr.aa + nr.bb + 1] <- ES.new * 
                1.01
        }
        if (increase == -1) {
            CES <- 1/ES.new - 1
            ans.all$ub[nr.var + nr.aa + nr.bb + 1] <- (1/ES.new)/1.01
        }
        ans.all$CES <- CES
        ans.all$pooled.sd <- pooled.sd
        return(ans.all)
    })
}



f.var.pooled <- function(ans.all) {
    if (exists("track")) 
        print("f.var.pooled")
    f.var <- function(dtype, xx, yy, sd2.log, nn) {
        if (dtype %in% c(10, 15, 250)) {
            SS <- sum(sd2.log * (nn - 1))
            df <- sum(nn - 1)
            var.within <- SS/df
        }
        else {
            if (dtype %in% c(1, 5)) 
                yy <- log(yy)
            if (dtype %in% c(26)) 
                yy <- sqrt(yy)
            mn <- tapply(yy, xx, mean)
            nn <- tapply(yy, xx, length)
            Vmn <- rep(mn[1], nn[1])
            if (dtype %in% c(1, 5)) 
                for (ii in 2:length(mn)) Vmn <- c(Vmn, rep(mn[ii], 
                  nn[ii]))
            resid <- yy - Vmn
            df <- sum(nn - 1)
            if (df == 0) 
                cat("\nconf intervals cannot be calculated, df = 0\n")
            var.within <- (sum(nn) - 1) * var(resid)/df
        }
        out.lst <- list(var.within = var.within, df = df)
        if (exists("track")) 
            print("f.var:   END")
        return(out.lst)
    }
    with(ans.all, {
        Vvar.within <- numeric()
        Vdf <- numeric()
        for (jj in 1:nr.bb) for (ii in 1:nr.aa) {
            x.part <- x[fct1 == ii & fct2 == jj]
            y.tmp <- yy[fct1 == ii & fct2 == jj]
            if (length(y.tmp) > 0) {
                if (dtype %in% c(10, 15, 250, 260)) {
                  sd2.log.part <- sd2.log[fct1 == ii & fct2 == 
                    jj]
                  nn.part <- nn[fct1 == ii & fct2 == jj]
                }
                else {
                  sd2.log.part <- NA
                  nn.part <- NA
                }
                out.lst <- f.var(dtype, xx = x.part, yy = y.tmp, 
                  sd2.log = sd2.log.part, nn.part)
                Vvar.within <- c(Vvar.within, out.lst$var.within)
                Vdf <- c(Vdf, out.lst$df)
            }
        }
        var.within.pooled <- (sum(Vdf * Vvar.within))/sum(Vdf)
        return(var.within.pooled)
    })
}



f.check.cc <- function(ans.all) {
    if (exists("track")) 
        print("f.check.cc")
    with(ans.all, {
        cc.OK <- T
        cc <- MLE[nr.var + nr.aa + nr.bb + 1]
        if (model.ans %in% c(4, 5, 9, 10, 14, 15, 24, 25, 52, 
            54)) {
            if (increase == 1) 
                if (cc < 0.02 + (1 + CES)) {
                  cat("\n\nthe value of parameter c is too close to CES\n                  this indicates (nonrandom) errors in the data,\n                  or too high value of CES")
                  cc.OK <- F
                  return(cc.OK)
                }
            if (increase == -1) {
                if (cc > -0.02 + (1 - abs(CES))) {
                  cat("\n\nthe value of parameter c is too close to CES\n                  this indicates (nonrandom) errors in the data,\n                  or too high value of CES")
                  cc.OK <- F
                }
                return(cc.OK)
            }
        }
        if (exists("track")) 
            print("f.check.cc    END")
        return(cc.OK)
    })
}



f.lines.con <- function(ans.all, display_plots = TRUE) {
    if (exists("track")) 
        print("f.lines.con")
    with(ans.all, {
        lines.plt.lst <- f.lines.plt.con(ans.all)
        nr.gr <- length(lines.plt.lst)
        for (ii in 1:nr.gr) {
            dfr <- lines.plt.lst[[ii]]
            x.tmp <- dfr$xline
            y.tmp <- dfr$expect
            if (display_plots) {
                lines(x.tmp[1:2], y.tmp[1:2], col = color[ii], pch = mark[ii], 
                    lty = 2)
                lines(x.tmp[-(1:2)], y.tmp[-(1:2)], col = color[ii], 
                    pch = mark[ii], lty = l.ty)
            }
        }
        ans.all$lines.plt.lst <- lines.plt.lst
        ans.all$nr.gr <- nr.gr
        if (exists("track")) 
            print("f.lines.con : END")
        return(ans.all)
    })
}



f.lines.plt.con <- function(ans.all) {
    if (exists("track")) 
        print("f.lines.plt.con")
    with(ans.all, {
        if (sum(is.na(regr.par)) > 0) {
            cat("\nregr.par contains NAs\n")
            print(regr.par)
            return(invisible())
        }
        if (plot.type == 12) 
            return(invisible())
        f.lines.tmp <- function(ans.all.tmp) {
            dfr <- with(ans.all.tmp, {
                nbins <- 1000
                nbins.0 <- nbins/10
                if (plot.type %in% c(2, 4, 8)) 
                  xline <- 10^(seq(from = log10(dum.contr), to = log10(max.x), 
                    length = nbins))
                else if (plot.type %in% 5:6) 
                  xline <- (seq(from = sqrt(min.x), to = sqrt(max.x), 
                    length = nbins))^2
                else xline <- seq(from = min.x, to = max.x, length = nbins)
                twice <- F
                x1.line <- NA
                x2.line <- NA
                if (model.ans == 42) {
                  if (m42.subgr == 1) {
                    x1.line <- xline
                    x2.line <- rep(0, length(xline))
                  }
                  if (m42.subgr == 2) {
                    x1.line <- rep(0, length(xline))
                    x2.line <- xline
                  }
                  if (m42.subgr == 3) {
                    x2.line <- xline
                    x1.line <- xline
                  }
                }
                if (incr.decr.no > 0) 
                  CES <- CES[1]
                expect <- f.expect.con(model.ans, xline, regr.par, 
                  fct1 = rep(1, length(xline)), fct2 = rep(1, 
                    length(xline)), fct3 = 1, CES = CES, twice = twice, 
                  trace = F, increase = increase, x.mn = x.mn, 
                  ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd, 
                  x1 = x1.line, x2 = x2.line, opposing = opposing, 
                  cc.inf = cc.inf)
                if (plot.type %in% c(2, 4, 8)) {
                  if (min.x == 0) {
                    xline <- log10(c(dum.contr, xline[nbins.0:nbins]))
                    if (model.ans == 42) {
                      x1.line <- log10(c(dum.contr, x1.line[nbins.0:nbins]))
                      x2.line <- log10(c(dum.contr, x2.line[nbins.0:nbins]))
                    }
                    expectdum <- f.expect.con(model.ans, 0, regr.par, 
                      fct1 = 1, fct2 = 1, fct3 = 1, CES = CES, 
                      twice = twice, increase = increase, x.mn = x.mn, 
                      ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd, 
                      x1 = 0, x2 = 0, opposing = opposing, cc.inf = cc.inf)
                    expect <- c(expectdum, expect[nbins.0:nbins])
                  }
                  else xline <- log10(xline)
                }
                if (plot.type %in% 5:6) 
                  xline <- sqrt(xline)
                if (plot.type %in% c(3, 4, 6)) 
                  expect <- log10(expect)
                if (plot.type %in% 7:8) 
                  expect <- sqrt(expect)
                dfr <- data.frame(xline = xline, expect = expect)
                if (plot.type %in% c(2, 4, 6) && min(x) == 0) 
                  attr(dfr, "dummy.index") <- 1
                else attr(dfr, "dummy.index") <- 0
                return(dfr)
            })
        }
        ans.all.tmp <- ans.all
        ans.all.tmp$dum.contr <- xy.lim[1]
        ans.all.tmp$min.x <- xy.lim[2]
        ans.all.tmp$max.x <- xy.lim[3]
        if (model.ans == 46) {
            ans.all.tmp$model.ans <- 15
            CED <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
            CED.ref <- CED[ref.lev]
            CED <- CED.ref/CED
            CED[ref.lev] <- CED.ref
            regr.par[(nr.aa + 1):(nr.aa + nr.bb)] <- CED
            ans.all.tmp$regr.par <- regr.par
        }
        if (sep) {
            regr.par.matr <- matrix(regr.par, nrow = 1)
        }
        else {
            if (model.ans %in% c(31, 33, 57)) 
                ans.all.tmp <- f.pars.m31(ans.all.tmp)
            else ans.all.tmp <- f.pars(ans.all.tmp)
            regr.par.matr <- ans.all.tmp$regr.par.matr
        }
        lines.plt.lst <- list()
        if ((model.ans != 6) && (nr.aa == 1) && (nr.bb == 1) && 
            (nr.cc == 1) && (nr.dd == 1) && (nr.var > 1)) {
            ans.all.tmp$regr.par <- regr.par.matr[1, ]
            dfr <- f.lines.tmp(ans.all.tmp)
            lines.plt.lst$all <- dfr
        }
        else {
            if (model.ans != 42) {
                nr.gr <- length(regr.par.matr[, 1])
                for (kk in 1:nr.gr) {
                  ans.all.tmp$regr.par <- regr.par.matr[kk, ]
                  if (ans.m6.sd == 2) 
                    ans.all.tmp$sign.q <- sign.q[kk]
                  dfr <- f.lines.tmp(ans.all.tmp)
                  if (nr.gr == 1) 
                    lines.plt.lst$all <- dfr
                  else lines.plt.lst[[gr.txt[kk]]] <- dfr
                }
            }
            if (model.ans == 42) {
                ans.all.tmp$m42.subgr <- 1
                dfr <- f.lines.tmp(ans.all.tmp)
                lines.plt.lst[[gr.txt[1]]] <- dfr
                ans.all.tmp$m42.subgr <- 2
                dfr <- f.lines.tmp(ans.all.tmp)
                lines.plt.lst[[gr.txt[2]]] <- dfr
            }
        }
        ans.all$nr.gr <- nr.gr
        if (exists("track")) 
            print("f.lines.plt.con:   END")
        return(lines.plt.lst)
    })
}



f.cedlines.con <- function(ans.all, display_plots = TRUE) {
    if (exists("track")) 
        print("f.cedlines.con")
    ans.all$cedes <- f.cedes.plt.con(ans.all)
    with(ans.all, {
        ES.x <- cedes$ES.x
        ES.y <- cedes$ES.y
        CED.x <- cedes$CED.x
        CED.y <- cedes$CED.y
        if (is.null(ans.all$regr.par.matr)) 
            regr.par.matr <- f.pars(ans.all)$regr.par.matr
        if (!is.matrix(regr.par.matr)) 
            regr.par.matr <- t(as.matrix(regr.par.matr))
        nr.gr <- length(regr.par.matr[, 1])
        if (length(xans) > 1) 
            nr.gr <- 1
        for (ii in (1:nr.gr)) {
            if (display_plots == TRUE) {
                lines(ES.x[ii, ], ES.y[ii, ], lty = 2)
                lines(CED.x[ii, ], CED.y[ii, ], lty = 2)
            }
        }
        if (exists("track")) 
            print("f.cedlines.con:  END")
        return(ans.all)
    })
}



f.cedes.plt.con <- function(ans.all) {
    if (exists("track")) 
        print("f.cedes.plt.con")
    if (is.null(ans.all$regr.par.matr)) 
        ans.all$regr.par.matr <- f.pars(ans.all)$regr.par.matr
    with(ans.all, {
        if (length(xans) > 1) 
            regr.par.matr <- f.pars(ans.all)$regr.par.matr
        if (model.ans == 46 && length(regr.par.matr[, 1]) > 1) {
            CED.ref <- CED[ref.lev]
            CED <- CED.ref/CED
            CED[ref.lev] <- CED.ref
            regr.par[(nr.aa + 1):(nr.aa + nr.bb)] <- CED
            ans.all$model.ans <- 15
            model.ans <- 15
            ans.all$regr.par <- regr.par
            regr.par.matr <- f.pars(ans.all)$regr.par.matr
        }
        if (!is.matrix(regr.par.matr)) 
            regr.par.matr <- t(as.matrix(regr.par.matr))
        nr.subgr <- length(regr.par.matr[, 1])
        if (ans.m6.sd == 2) {
            regr.par <- regr.par.matr[1, ]
            GCES <- CES
            qq <- regr.par[4]
            if (model.ans == 6) {
                if (is.na(CED)) 
                  cat("\n\nATTENTION: CED not yet calculated ......  \n\n")
                CED.0 <- CED
            }
            if (model.ans == 47) 
                CED.0 <- regr.par[3]
        }
        if (nr.bb == 1) 
            CED <- rep(CED, nr.aa)
        dum.contr <- xy.lim[1]
        low.y <- xy.lim[4]
        if (is.na(low.y)) 
            low.y <- min(y)
        if (low.y == 9999) 
            low.y <- min(y)
        ES.x <- numeric()
        ES.y <- numeric()
        CED.x <- numeric()
        CED.y <- numeric()
        for (jj in 1:nr.subgr) {
            if (ans.m6.sd == 2) {
                regr.par <- regr.par.matr[jj, ]
                sign.q.0 <- sign.q[jj]
                CES.0 <- exp(qq * sqrt(regr.par[1]))^GCES - 1
                ES.0 <- exp(log(regr.par[2]) + sign.q.0 * log(CES.0 + 
                  1))
            }
            else {
                CED.0 <- CED[jj]
                CES.0 <- CES
                if (incr.decr.no > 0) {
                  CES.0 <- abs(CES[1])
                  if (regr.par[3] > 1) 
                    CES.0 <- -CES.0
                }
            }
            ES.ced <- f.expect.con(model.ans, CED.0, regr.par.matr[jj, 
                ], fct1 = 1, fct2 = 1, fct5 = 1, CES = CES.0, 
                increase = increase, ref.lev = ref.lev, cc.inf = cc.inf)
            ES.backgr <- f.expect.con(model.ans, 0, regr.par.matr[jj, 
                ], fct1 = 1, fct2 = 1, fct5 = 1, CES = CES.0, 
                increase = increase, ref.lev = ref.lev, cc.inf = cc.inf)
            ES.0 <- ES.backgr * (CES.0 + 1)
            if (incr.decr.no == 0) 
                if (model.ans != 49) 
                  if (abs(ES.ced - ES.backgr * (CES + 1)) > 0.3) {
                    cat("\n response at CED is", ES.ced)
                    cat("\n response at CES is", ES.0)
                    cat("\n\nATTENTION: value of BMD appears incorrect\n")
                    cat("this indicates a problem in your analysis\n\n")
                    cat("in subgroup", jj, "\n\n")
                  }
            if (plot.type %in% c(1, 2, 5)) 
                ES.y <- rbind(ES.y, rep(ES.0, 2))
            if (plot.type %in% c(3, 4, 6)) 
                ES.y <- rbind(ES.y, rep(log10(ES.0), 2))
            if (plot.type %in% 7:8) 
                ES.y <- rbind(ES.y, rep(sqrt(ES.0), 2))
            if (plot.type %in% c(1, 3, 7)) 
                ES.x <- rbind(ES.x, c(min(x), CED.0))
            if (plot.type %in% c(2, 4, 8)) 
                ES.x <- rbind(ES.x, c(log10(dum.contr), log10(CED.0)))
            if (plot.type %in% c(5, 6)) 
                ES.x <- rbind(ES.x, c(0, sqrt(CED.0)))
            if (plot.type %in% c(1, 2, 5)) 
                CED.y <- rbind(CED.y, c(low.y, ES.0))
            if (plot.type %in% c(3, 4, 6)) 
                CED.y <- rbind(CED.y, c(log10(low.y), log10(ES.0)))
            if (plot.type %in% c(7, 8)) 
                CED.y <- rbind(CED.y, c(log10(low.y), sqrt(ES.0)))
            if (plot.type %in% c(1, 3, 7)) 
                CED.x <- rbind(CED.x, rep(CED.0, 2))
            if (plot.type %in% c(2, 4, 8)) 
                CED.x <- rbind(CED.x, rep(log10(CED.0), 2))
            if (plot.type %in% c(5, 6)) 
                CED.x <- rbind(CED.x, rep(sqrt(CED.0), 2))
        }
        ES.x <- matrix(ES.x, ncol = 2)
        ES.y <- matrix(ES.y, ncol = 2)
        CED.x <- matrix(CED.x, ncol = 2)
        CED.y <- matrix(CED.y, ncol = 2)
        if (0) 
            if (plot.type %in% c(3, 4, 6) && !is.na(xy.lim[6])) {
                CED.y <- CED.y + xy.lim[6]
                ES.y <- ES.y + xy.lim[6]
            }
        cedes.lst <- list()
        cedes.lst$CED.x = CED.x
        cedes.lst$CED.y = CED.y
        cedes.lst$ES.x = ES.x
        cedes.lst$ES.y = ES.y
        if (0) {
            print("f.cedes.plt.con")
            print(CED.x)
            print(CED.y)
            print(ES.x)
            print(ES.y)
        }
        if (exists("track")) 
            print("f.cedes.plt.con: END")
        return(cedes.lst)
    })
}



f.show.con <- function(ans.all) {
    if (exists("track")) 
        print("f.show.con")
    with(ans.all, {
        nl <- "\n"
        tb <- "  "
        sp <- ": "
        if (length(ans.all$fct4) > 0) 
            nr.cc <- max(fct4)
        if (model.ans %in% c(3, 8, 13, 18, 23, 51, 53)) 
            nr.cc <- 0
        if (length(ans.all$fct5) > 0) 
            nr.dd <- max(fct5)
        if (fit.ans == 1) {
            loglik.txt <- "loglik"
        }
        if (fit.ans == 2) {
            loglik.txt <- "SP"
            loglik <- SP
            nr.var <- 0
        }
        if (length(MLE) < 22) {
            if (dtype %in% c(5, 15)) {
                MLE.tmp <- c(inter.var, intra.var, MLE[-1])
                text.par <- c("var.inter", "var.intra", text.par[-1])
            }
            else MLE.tmp <- MLE
            results <- c(round(loglik, 2), round(aic, 2), signif(MLE.tmp, 
                4))
            txt <- c(loglik.txt, "AIC", text.par)
        }
        else {
            cat("\nATTENTION: Number of parameters too large for legend in plot\n")
            if (!model.ans %in% c(31, 33, 57)) {
                results <- c(round(loglik, 2), nr.var, nr.aa, 
                  nr.bb, nr.cc, nr.dd, signif(MLE[c(1, nr.var + 
                    1, nr.var + nr.aa + 1, nr.var + nr.aa + nr.bb + 
                    1, nr.var + nr.aa + nr.bb + nr.cc + 1)], 
                    3))
                txt <- c(loglik.txt, "nr.var", "nr.aa", "nr.bb", 
                  "nr.cc", "nr.dd", text.par[c(1, nr.var + 1, 
                    nr.var + nr.aa + 1, nr.var + nr.aa + nr.bb + 
                      1, nr.var + nr.aa + nr.bb + nr.cc + 1)])
            }
            if (model.ans %in% c(31, 33, 57)) {
                nr.b1 <- nr.dd
                nr.b2 <- nr.bb
                nr.dd <- 1
                ra.var <- 1
                ra.aa <- nr.var + 1
                ra.b1 <- nr.var + nr.aa + 1
                ra.b2 <- nr.var + nr.aa + nr.b1 + 1
                ra.c1 <- nr.var + nr.aa + nr.b1 + nr.b2 + 1
                ra.c2 <- nr.var + nr.aa + nr.b1 + nr.b2 + nr.cc + 
                  1
                ra.dd <- nr.var + nr.aa + nr.b1 + nr.b2 + nr.cc + 
                  2
                if (model.ans == 33) 
                  ra.dd <- c(ra.dd, nr.var + nr.aa + nr.b1 + 
                    nr.b2 + nr.cc + 3)
                if (model.ans == 57) {
                  ra.dd <- c(ra.c2, ra.dd)
                  ra.c2 <- NULL
                }
                mle.tmp <- MLE[c(ra.var, ra.aa, ra.b1, ra.b2, 
                  ra.c1, ra.c2, ra.dd)]
                text.tmp <- text.par[c(ra.var, ra.aa, ra.b1, 
                  ra.b2, ra.c1, ra.c2, ra.dd)]
                results <- c(round(loglik, 2), nr.var, nr.aa, 
                  nr.b1, nr.b2, nr.cc, 1, nr.dd, signif(mle.tmp), 
                  3)
                txt <- c(loglik.txt, "nr.var", "nr.aa", "nr.b1", 
                  "nr.b2", "nr.c1", "nr.c2", "nr.dd", text.tmp)
            }
            if (model.ans %in% c(13, 23, 51, 53)) {
                results <- results[-c(5, 10)]
                txt <- txt[-c(5, 11)]
            }
            if (length(xans) > 1) {
                RPF <- MLE[(length(MLE) - nr.dosecol + 2):length(MLE)]
                RPF.txt <- text.par[(length(text.par) - nr.dosecol + 
                  2):length(text.par)]
                results <- c(results, signif(RPF, 3))
                txt <- c(txt, RPF.txt)
            }
            results <- results[!is.na(results)]
            txt <- txt[!is.na(txt)]
        }
        show <- paste(" version:", PRversion[1])
        for (i in 1:length(txt)) show <- paste(show, nl, txt[i], 
            tb, results[i])
        if (!is.na(CED[1]) || CED.model) {
            if (nr.bb == 1) 
                fct2.txt <- NULL
            CES.txt <- "CES"
            if (NES.ans == 2 || model.ans == 47) 
                CES.txt <- "adjusted CES"
            CED.txt <- "CED"
            if (length(ans.all$m6.ED50) != 0) 
                if (m6.ED50) 
                  CED.txt <- "ED50"
                else m6.ED50 <- FALSE
            if (incr.decr.no == 0) {
                CES.tmp <- round(CES[1], 3)
                if (CES < 0.01) 
                  CES.tmp <- signif(CES, 3)
                if ((CES != 0 || m6.ED50) && !is.na(CED[1])) 
                  show <- paste(show, nl, CES.txt, tb, CES.tmp)
                if (!CED.model) {
                  CED.uniq <- unique(CED)
                  if (model.ans == 43) {
                    nr.bb <- length(CED.uniq)
                    fct2.txt <- fct1.txt
                  }
                  for (ii in (1:nr.bb)) show <- paste(show, nl, 
                    paste(CED.txt, fct2.txt[ii], sep = "-"), 
                    signif(CED.uniq[ii], 4))
                }
                else {
                  if (model.ans == 47) 
                    show <- paste(show, nl, "adjusted CES", tb, 
                      CES.tmp)
                }
            }
            else CES <- CES[1]
            if (CES == 0) 
                show <- paste(show, nl, "BMR", BMR, "SD")
            if (length(conf.int) > 0) 
                if (any(!is.na(conf.int))) {
                  if (boot) 
                    bt <- "bt"
                  else bt <- ""
                  BMDL.txt <- paste("CEDL", bt, sep = "")
                  BMDU.txt <- paste("CEDU", bt, sep = "")
                  if (model.ans == 46) {
                    nr.bb <- length(fct2.txt)
                    BMDL.txt <- rep("RPFL", length(fct2.txt))
                    BMDU.txt <- rep("RPFU", length(fct2.txt))
                    BMDL.txt[ref.lev] <- "CEDL"
                    BMDU.txt[ref.lev] <- "CEDU"
                    BMDL.txt <- paste(BMDL.txt, bt, sep = "")
                    BMDU.txt <- paste(BMDU.txt, bt, sep = "")
                  }
                  if (model.ans == 46) {
                    for (jj in (1:nr.bb)) {
                      show <- paste(show, nl, paste(BMDL.txt[jj], 
                        fct2.txt[jj], sep = "-"), tb, signif(conf.int[jj, 
                        1], 4))
                      show <- paste(show, nl, paste(BMDU.txt[jj], 
                        fct2.txt[jj], sep = "-"), tb, signif(conf.int[jj, 
                        2], 4))
                    }
                  }
                  else if ((nr.aa != nr.bb) && (nr.bb > 1) && 
                    (nr.aa > 1)) 
                    for (ii in (1:nr.bb)) {
                      show <- paste(show, nl, paste(BMDL.txt, 
                        fct2.txt[ii], sep = "-"), tb, signif(conf.int[ii, 
                        1], 4))
                      show <- paste(show, nl, paste(BMDU.txt, 
                        fct2.txt[ii], sep = "-"), tb, signif(conf.int[ii, 
                        2], 4))
                    }
                  else if ((model.ans == 16) | (nr.bb == 1)) {
                    show <- paste(show, nl, BMDL.txt, tb, signif(conf.int[1, 
                      1], 4))
                    show <- paste(show, nl, BMDU.txt, tb, signif(conf.int[1, 
                      2], 4))
                  }
                  else if ((model.ans != 16) & (nr.aa > 1)) {
                    for (ii in (1:nr.aa)) {
                      show <- paste(show, nl, paste(BMDL.txt, 
                        fct2.txt[ii], sep = "-"), tb, signif(conf.int[ii, 
                        1], 4))
                      show <- paste(show, nl, paste(BMDU.txt, 
                        fct2.txt[ii], sep = "-"), tb, signif(conf.int[ii, 
                        2], 4))
                    }
                  }
                  else for (ii in (1:nr.bb)) {
                    show <- paste(show, nl, paste(BMDL.txt, fct2.txt[ii], 
                      sep = "-"), tb, signif(conf.int[ii, 1], 
                      4))
                    show <- paste(show, nl, paste(BMDU.txt, fct2.txt[ii], 
                      sep = "-"), tb, signif(conf.int[ii, 2], 
                      4))
                  }
                  if (nruns > 0) 
                    show <- paste(show, nl, paste("nruns:", nruns))
                }
        }
        if (incr.decr.no == 0) 
            if (model.ans %in% c(12:15, 22:25, 46, 51, 52, 53, 
                54)) {
                CED.tmp <- MLE[(nr.var + nr.aa + 1):(nr.var + 
                  nr.aa + nr.bb)]
                cc0 <- NA
                dd0 <- NA
                if (model.ans %in% c(13, 23, 51, 53)) {
                  if (increase == 1) 
                    cc0 <- cc.inf
                  if (increase == -1) 
                    cc0 <- 1/cc.inf
                  show <- paste(show, nl, paste("c-  fixed at", 
                    cc0))
                  dd0 <- MLE[(nr.var + nr.aa + nr.bb + 1):(nr.var + 
                    nr.aa + nr.bb + nr.dd)]
                }
                if (model.ans %in% c(14:15, 24:25, 46, 52, 54)) 
                  cc0 <- MLE[(nr.var + nr.aa + nr.bb + 1):(nr.var + 
                    nr.aa + nr.bb + nr.cc)]
                if (model.ans %in% c(15, 25, 46, 52, 54)) 
                  dd0 <- MLE[(nr.var + nr.aa + nr.bb + nr.cc + 
                    1):(nr.var + nr.aa + nr.bb + nr.cc + nr.dd)]
                bb.tmp <- f.bb.con(model.ans, cc0, dd0, CED.tmp, 
                  CES, CED.ref = CED.tmp[ref.lev], ref.lev = ref.lev, 
                  cont = TRUE)
                if (length(MLE) < 10) 
                  for (i in 1:length(bb.tmp)) show <- paste(show, 
                    "\n b: ", signif(bb.tmp[i], 4))
            }
        if (dtype %in% c(5, 15)) {
            show <- paste(show, nl, "inter.var:", signif(inter.var, 
                2))
            show <- paste(show, nl, "intra.var:", signif(intra.var, 
                2))
        }
        show <- paste(show, nl, "conv", sp, converged)
        show <- paste(show, nl, "scaling factor on x", sp, sf.x)
        if (any(Vdetlim > 0)) {
            if (length(unique(Vdetlim)) == 1) 
                show <- paste(show, nl, "detlim", sp, Vdetlim[1])
            else show <- paste(show, nl, "detlim", sp, min(Vdetlim, 
                na.rm = T), "-", max(Vdetlim, na.rm = T))
        }
        if (!is.na(cens.up)) 
            show <- paste(show, nl, "cens.up", sp, cens.up)
        show <- paste(show, nl, "dtype", sp, dtype)
        if (!gui) {
            sel.tmp <- character(0)
            for (ii in 1:length(select.name)) sel.tmp <- paste(sel.tmp, 
                "selected", sp, select.name[ii])
            if (length(sel.tmp) > 0) 
                show <- paste(show, nl, sel.tmp)
            else show <- paste(show, nl, "selected: all")
            rem.tmp <- character()
            if (is.na(remove.name[1])) 
                show <- paste(show, nl, "removed: NA")
            else if (remove.name[1] == "none") 
                show <- paste(show, nl, "removed: none")
            else {
                for (ii in 1:length(remove.name)) rem.tmp <- paste(rem.tmp, 
                  " ", remove.name[ii])
                if (length(rem.tmp) >= 1) {
                  show <- paste(show, nl, "removed:")
                  show <- paste(show, rem.tmp)
                }
            }
        }
        if (quick.ans == 1) 
            show <- paste(show, nl, factor.name)
        if (length(covar.txt) > 1) 
            show <- paste(show, nl, paste("covariate:", covar.name, 
                "\n"))
        if (cond.ans != 3) 
            show <- paste(show, "\n fit cond:", cond.ans)
        if (exists("track")) 
            print("f.show.con:  END")
        return(show)
    })
}



f.plot.all <- function(ans.all, sep = F, bootstrap = F, new.window = TRUE, 
    no.show = FALSE, display_plots = TRUE) {
    if (exists("track")) 
        print("f.plot.all")
    if (ans.all$plot.type == 0) 
        return(ans.all)
    if (ans.all$tans != 0) {
        ans.all$y.lim <- f.plot.cxt(ans.all, sep = sep)
        return(ans.all)
    }
    if (is.na(ans.all$xy.lim[1]) || ans.all$xy.lim[1] == 0) {
        x <- ans.all$x
        dum.contr <- min(x[x > 0])/5
        ans.all$xy.lim[1] <- dum.contr
    }
    if (ans.all$dtype == 3 && !ans.all$combi.ans && (ans.all$plot.type %in% 
        5:6)) {
        if (ans.all$quick.ans == 1) 
            ans.all$categ.ans <- menu(c("all categories in one plot", 
                "separate plots for each category"), title = paste("\n\nHow do you want to plot the categories?\n"))
        else ans.all$categ <- 1
        if (ans.all$nr.gr > 1) 
            sep <- TRUE
        if (ans.all$displ.no > 0) {
            sep <- TRUE
            ans.all$fct2 <- ans.all$displ.fact
        }
        if (ans.all$categ.ans == 2) 
            ans.all$cex.1 <- 2
        ans.all$sep <- sep
    }
    if (length(ans.all$xans) > 1 && ans.all$nr.aa > 1) {
        sep <- TRUE
        cat("\n\nATTENTION: plots are provided for each subgroup, as global plot is not correct \n for dose additon model with covariate on parameter a\n\n")
        f.press.key.to.continue()
    }
    if (sep) {
        f.plot.sep(ans.all)
        return(ans.all)
    }
    with(ans.all, {
        if (cont) {
            if (new.window)
                if (display_plots) {
                  f.graph.window(1)
                }
            if (plot.type %in% 10:11) {
                if (is.na(regr.resid[1])) 
                  cat("\nno residuals available \n")
                else ans.all$yy <- regr.resid
                ans.all$CI.plt <- FALSE
            }
            ans.all$heading <- ""
            ans.all <- f.plot.con(ans.all, display_plots = display_plots)
            if (model.ans != 11 && plot.type != 9) 
                ans.all <- f.lines.con(ans.all, display_plots = display_plots)
            if (plot.type < 9 && !is.na(CED[1]) && model.ans != 
                16 && fitted) 
                ans.all <- f.cedlines.con(ans.all, display_plots = display_plots)
            if (is.na(ans.all$MLE[1])) 
                ans.all$MLE <- ans.all$par.start
            if (no.show) 
                ans.all$show <- ""
            else ans.all$show <- f.show.con(ans.all)
            f.mtext(ans.all, display_plots = display_plots)
        }
        if (!cont) {
            if (!combi.ans) {
                if (!bootstrap && new.window) 
                  if (display_plots) {
                    f.graph.window(1)
                  }
                if (dtype %in% c(4, 6, 84) && model.type == 2 && 
                  plot.type %in% c(7:10)) {
                  ans.all.plt <- ans.all
                  out.list <- f.lines.cat(ans.all.plt)
                  ans.all.plt$y.lim <- out.list$y.lim
                  return(ans.all)
                }
                else if (dtype %in% c(4, 6, 84) || (dtype == 
                  2 && plot.type < 3)) {
                  ans.all <- f.plot.frq(ans.all)
                  if (length(unique(fct3)) > 1) 
                    ans.all$plotcurve <- TRUE
                  if (!(model.ans == 14 && (model.type == 1))) {
                    f.lines.frq(ans.all)
                  }
                  if (!all(is.na(ans.all$CED.matr))) 
                    f.cedlines.bin(ans.all)
                }
                else {
                  ans.all.plt <- ans.all
                  out.list <- f.lines.cat(ans.all.plt)
                  ans.all$y.lim.plt <- out.list$y.lim
                }
                if (ans.all$categ.ans < 2) 
                  f.mtext(ans.all)
            }
            if (combi.ans) {
                if (display_plots) {
                  f.graph.window(4)
                }
                ans.all.plt <- ans.all
                ans.all.plt$modelname <- ""
                ans.all.plt$main <- ""
                ans.all.plt$plot.type <- 4 + ans.scale
                if (dtype %in% c(2, 3, 4, 6)) 
                  ans.all.plt$ced.lines <- TRUE
                out.list <- f.lines.cat(ans.all.plt)
                ans.all.plt$y.lim <- out.list$y.lim
                if (no.show) 
                  ans.all$show <- ""
                else ans.all.plt$show <- f.show.cat(ans.all.plt)
                title(modelname, col = color[1], cex.main = 0.9, 
                  font.main = 1)
                f.mtext(ans.all.plt)
                ans.all.plt$plot.type <- 5 + ans.scale
                out.list <- f.lines.cat(ans.all.plt)
                ans.all.plt$y.lim.plt <- out.list$y.lim
                title(modelname, col = color[1], cex.main = 0.9, 
                  font.main = 1)
                if (0) {
                  ans.all.plt$plot.type <- 9 + ans.scale
                  ans.all.plt$modelname <- ""
                  ans.all.plt$ced.lines <- F
                  out.list <- f.lines.cat(ans.all.plt)
                  ans.all.plt$y.lim <- out.list$y.lim
                  if (!is.na(CED.matr[1])) {
                    th.par.cum <- cumsum(th.par)
                    CED.vect <- CED.matr[1, ]
                    if (ans.scale == 1) {
                      CED.vect <- log10(CED.vect)
                      x <- log10(x[x > 0]/sf.x)
                    }
                    for (ii in 1:length(CED.vect)) lines(c(min(x), 
                      CED.vect[ii]), rep(th.par.cum[ii], 2), 
                      lty = 2)
                    y.leg.2 <- paste(y.leg, " (latent var.uu)")
                    if (dtype == 3 & nr.gr == 1) {
                      CED.tmp <- round(CED.vect, rep(4, length(CED.vect)))
                      CED.tmp <- CED.tmp[CED.tmp < max(x)]
                      CED.left <- c(CED.tmp[1]/100, CED.tmp)
                      CED.right <- c(CED.tmp, max(x))
                      CED.mid <- (CED.left + CED.right)/2
                      if (ans.scale == 1) 
                        CED.mid <- 10^(CED.mid)
                      uu <- f.expect.con(model.ans, CED.mid, 
                        regr.par, fct1 = 1, fct2 = 1, CES = CES, 
                        increase = increase, cc.inf = cc.inf)
                      zz <- logb(uu)
                      if (ans.scale == 1) 
                        CED.mid <- log10(CED.mid)
                      for (ii in 1:length(scores.orig)) text(CED.mid[ii], 
                        zz[ii], scores.orig[ii])
                    }
                  }
                }
                th.par.cum <- cumsum(th.par)
                if (mode(CED.matr) == "NULL") 
                  CED.vect <- 0
                else if (is.na(CED.matr[1])) 
                  CED.vect <- 0
                else CED.vect <- CED.matr[1, ]
                if (dtype %in% 2:3) {
                  ans.all.plt$plot.type <- 5 + ans.scale
                }
            }
        }
        if (0) {
            if (model.ans == 46 && dtype == 3) 
                ans.all <- ans.all.saved
        }
        if (exists("track")) 
            print("f.plot.all : END")
        return(ans.all)
    })
}



f.move.sublist <- function(lst, sub) {
    if (exists("track")) 
        print("f.move.sublist")
    lst$model.ans <- sub$model.ans
    if (length(sub$model.type) > 0) 
        lst$model.type <- sub$model.type
    lst$modelname <- sub$modelname
    lst$MLE <- sub$MLE
    lst$MLE.all <- sub$MLE.all
    if (length(sub$CED) > 0) {
        lst$CED <- sub$CED
        lst$CED.matr <- sub$CED.matr
    }
    lst$response.matr <- sub$response.matr
    lst$report.pars <- sub$report.pars
    if (length(lst$xans) > 1) {
        lst$x <- sub$x
        lst$xy.lim <- sub$xy.lim
    }
    if (lst$NES.ans == 2) 
        lst$CES <- sub$CES
    lst$par.start <- sub$MLE
    lst$loglik.old <- sub$loglik
    lst$regr.par <- sub$regr.par
    lst$regr.par.matr <- sub$regr.par.matr
    if (!lst$cont) {
        lst$th.par <- sub$th.par
        lst$sig.par <- sub$sig.par
        lst$nth <- sub$nth
        lst$sens.lev <- sub$sens.lev
        lst$rank.low <- sub$rank.low
    }
    if (length(sub$nrp) > 0) 
        lst$nrp <- sub$nrp
    lst$npar <- sub$npar
    lst$npar.aic <- sub$npar.aic
    lst$text.par <- sub$text.par
    lst$lower <- sub$lower
    lst$upper <- sub$upper
    lst$lb <- sub$lb
    lst$ub <- sub$ub
    lst$loglik <- sub$loglik
    lst$aic <- sub$aic
    lst$converged <- sub$converged
    lst$conf.int <- sub$conf.int
    if (lst$dtype == 3) 
        lst$conf.int.matr <- sub$conf.int.matr
    lst$fct1 <- sub$fct1
    lst$fct2 <- sub$fct2
    lst$fct3 <- sub$fct3
    lst$fct4 <- sub$fct4
    lst$fct5 <- sub$fct5
    lst$nr.aa <- max(sub$fct1)
    lst$nr.bb <- max(sub$fct2)
    lst$nr.var <- max(sub$fct3)
    lst$nr.cc <- max(sub$fct4)
    lst$nr.dd <- max(sub$fct5)
    lst$nr.gr <- sub$nr.gr
    lst$show <- sub$show
    lst$trend <- sub$trend
    lst$regr.resid <- sub$regr.resid
    if (lst$dtype %in% c(5, 15)) {
        lst$inter.var <- sub$inter.var
        lst$intra.var <- sub$intra.var
    }
    lst$gr.txt <- sub$gr.txt
    lst$lines.plt.lst <- sub$lines.plt.lst
    lst$cedes.lst <- sub$cedes.lst
    if (lst$dtype %in% c(1, 5)) 
        lst$outliers <- sub$outliers
    if (exists("track")) 
        print("f.move.sublist:  END")
    return(lst)
}



f.ced.ma <- function(ans.all, xline, first.call, display_plots = TRUE) {
    if (exists("track")) 
        print("f.ced.ma")
    with(ans.all, {
        length.xx.interpol <- 1000
        ans.all$nr.models <- nr.models
        ans.all$model.list <- model.list
        nr.gr <- max(covariate)
        if (!cont) {
            Vaic <- Vaic[-(1:2)]
        }
        min.aic <- min(Vaic)
        Vaic.adj <- Vaic - min.aic
        weight.denom <- sum(exp(-0.5 * Vaic.adj))
        Vweight <- exp(-0.5 * Vaic.adj)/weight.denom
        Vweight.dfr = data.frame(model = model.list[1:nr.models], 
            weight = round(Vweight, 4))
        if (first.call) {
            xx.ma <- numeric()
            yy.ma <- numeric()
            nn.ma <- numeric()
            covar.ma <- numeric()
            if (dtype %in% c(1, 10, 25, 250, 26, 260)) 
                sigma.ma <- numeric()
            if (dtype %in% c(5, 15)) {
                inter.sigma.ma <- numeric()
                intra.sigma.ma <- numeric()
            }
            for (group in 1:nr.gr) {
                x.subgr <- x[covariate == group]
                if (dtype %in% c(1, 5, 25, 26)) 
                  x.subgr <- unique(x.subgr)
                y.subgr <- numeric(length(x.subgr))
                if (dtype %in% c(1, 10, 25, 250, 26, 260)) 
                  var.subgr <- numeric(1)
                if (dtype %in% c(5, 15)) {
                  inter.var.subgr <- numeric(1)
                  intra.var.subgr <- numeric(1)
                }
                if (dtype %in% c(1, 5, 25, 26)) {
                  nn.table <- table(x, covariate)
                  n.subgr <- nn.table[, group]
                  n.subgr <- n.subgr[n.subgr > 0]
                }
                else n.subgr <- nn[covariate == group]
                Vmodel.ans <- numeric()
                Vmodelname <- numeric()
                for (ii in 1:nr.models) {
                  model.ans <- ans.all[[model.list[ii]]]$model.ans
                  regr.par.matr <- ans.all[[model.list[ii]]]$regr.par.matr
                  nr.col <- length(regr.par.matr[1, ])
                  if (length(regr.par.matr[, 1]) == 1 && nr.gr > 
                    1) 
                    regr.par.matr <- matrix(regr.par.matr, ncol = nr.col, 
                      nrow = nr.gr, byrow = TRUE)
                  regr.par <- regr.par.matr[group, ]
                  Vmodel.ans <- c(Vmodel.ans, model.ans)
                  Vmodelname <- c(Vmodelname, ans.all[[model.list[ii]]]$modelname)
                  if (!cont) {
                    if (ii > (nr.models - 2)) {
                      th.par <- ans.all[[model.list[ii]]]$th.par
                      sig.par <- ans.all[[model.list[ii]]]$sig.par
                      y.expect <- f.expect.cat(model.type = 2, 
                        model.ans, x.subgr, regr.par, th.par, 
                        sig.par, dtype = dtype, CES = CES, ces.ans = ces.ans, 
                        fct1 = 1, fct2 = 1, cc.inf = cc.inf)
                    }
                    else {
                      y.expect <- f.expect.bin(model.ans, x.subgr, 
                        regr.par, dtype = dtype, CES = CES, ces.ans = ces.ans, 
                        fct1 = 1, fct2 = 1, output = output)
                    }
                  }
                  if (cont) 
                    y.expect <- f.expect.con(model.ans, x.subgr, 
                      regr.par, CES = CES, increase = increase, 
                      cc.inf = cc.inf)
                  if (all(!is.na(y.expect))) {
                    if (logit.ma) 
                      y.expect <- log(y.expect)/(1 - log(y.expect))
                    if (cont) 
                      y.expect <- log(y.expect)
                    y.subgr <- y.subgr + Vweight[ii] * y.expect
                    if (!cont && !logit.ma) 
                      y.subgr[y.subgr > 1] <- 1
                  }
                  else cat("\n\nATTENTION:   NAs in model average response \n\n")
                  if (dtype %in% c(1, 10, 25, 250, 26, 260)) {
                    MLE <- ans.all[[model.list[ii]]]$MLE
                    if (max(ans.all[[model.list[ii]]]$fct3) > 
                      1) 
                      res.var <- MLE[group]
                    else res.var <- MLE[1]
                    var.subgr <- var.subgr + Vweight[ii] * res.var
                  }
                  if (dtype %in% c(5, 15)) {
                    inter.var <- ans.all[[model.list[ii]]]$inter.var
                    intra.var <- ans.all[[model.list[ii]]]$intra.var
                    inter.var.subgr <- inter.var.subgr + Vweight[ii] * 
                      inter.var
                    intra.var.subgr <- intra.var.subgr + Vweight[ii] * 
                      intra.var
                  }
                }
                if (logit.ma) 
                  y.subgr <- exp(y.subgr/(y.subgr + 1))
                if (cont) 
                  y.subgr <- exp(y.subgr)
                yy.ma <- c(yy.ma, y.subgr)
                xx.ma <- c(xx.ma, x.subgr)
                nn.ma <- c(nn.ma, n.subgr)
                if (dtype == 2) 
                  nn.ma <- 1
                covar.ma <- c(covar.ma, rep(group, length(x.subgr)))
                if (dtype %in% c(1, 10, 25, 250, 26, 260)) 
                  sigma.ma <- c(sigma.ma, sqrt(var.subgr))
                if (dtype %in% c(5, 15)) {
                  inter.sigma.ma <- c(inter.sigma.ma, sqrt(inter.var.subgr))
                  intra.sigma.ma <- c(intra.sigma.ma, sqrt(intra.var.subgr))
                }
                if (exists("track")) 
                  print("f.ced.ma:   END 1")
            }
            ans.all$MA <- list(yy.ma = yy.ma, xx.ma = xx.ma, 
                nn.ma = nn.ma, covar.ma = covar.ma, Vweight = Vweight.dfr, 
                Vmodel.ans = Vmodel.ans, Vmodelname = Vmodelname)
            if (dtype %in% c(1, 10, 25, 26, 250, 260)) 
                ans.all$MA$sigma.ma <- sigma.ma
            if (dtype %in% c(5, 15)) {
                ans.all$MA$inter.sigma.ma <- inter.sigma.ma
                ans.all$MA$intra.sigma.ma <- intra.sigma.ma
            }
            return(ans.all$MA)
        }
        if (!first.call) {
            Vced.ma <- numeric()
            for (group in 1:nr.gr) {
                yy.ma <- numeric(length.xx.interpol + 1)
                yy.ma.plt <- numeric(length(xline))
                xx.interpol <- f.xx.ma(ans.all, group = group, 
                  length.xx = length.xx.interpol)
                for (ii in 1:nr.models) {
                  model.ans <- ans.all[[model.list[ii]]]$model.ans
                  regr.par.matr <- ans.all[[model.list[ii]]]$regr.par.matr
                  nr.col <- length(regr.par.matr[1, ])
                  if (length(regr.par.matr[, 1]) == 1 && nr.gr > 
                    1) 
                    regr.par.matr <- matrix(regr.par.matr, ncol = nr.col, 
                      nrow = nr.gr, byrow = TRUE)
                  regr.par <- regr.par.matr[group, ]
                  if (!cont) {
                    if (ii > (nr.models - 2)) {
                      th.par <- ans.all[[model.list[ii]]]$th.par
                      sig.par <- ans.all[[model.list[ii]]]$sig.par
                      y.expect <- f.expect.cat(model.type = 2, 
                        model.ans, xx.interpol, regr.par, th.par, 
                        sig.par, dtype = dtype, CES = CES, ces.ans = ces.ans, 
                        fct1 = 1, fct2 = 1, cc.inf = cc.inf)
                      y.expect.plt <- f.expect.cat(model.type = 2, 
                        model.ans, xline, regr.par, th.par, sig.par, 
                        dtype = dtype, CES = CES, ces.ans = ces.ans, 
                        fct1 = fct1, fct2 = 1, cc.inf = cc.inf)
                    }
                    else {
                      y.expect <- f.expect.bin(model.ans, xx.interpol, 
                        regr.par, dtype = dtype, CES = CES, ces.ans = ces.ans, 
                        fct1 = 1, fct2 = 1, output = output)
                      y.expect.plt <- f.expect.bin(model.ans, 
                        xline, regr.par, dtype = dtype, CES = CES, 
                        ces.ans = ces.ans, fct1 = 1, fct2 = 1, 
                        output = output)
                    }
                  }
                  if (cont) {
                    y.expect <- f.expect.con(model.ans, xx.interpol, 
                      regr.par, CES = CES, increase = increase, 
                      cc.inf = cc.inf)
                    y.expect.plt <- f.expect.con(model.ans, xline, 
                      regr.par, CES = CES, increase = increase, 
                      cc.inf = cc.inf)
                    y.expect <- log(y.expect)
                    y.expect.plt <- log(y.expect.plt)
                  }
                  if (0) {
                    print(regr.par)
                    print("weights")
                    print(Vweight)
                    print(cc.inf)
                    print(increase)
                    print(CES)
                    print("expect")
                    print(y.expect[c(1, 100, 500)])
                    print("xx")
                    print(xx.interpol[c(1, 100, 500)])
                    x.pr()
                  }
                  if (all(!is.na(y.expect))) {
                    yy.ma <- yy.ma + Vweight[ii] * y.expect
                    yy.ma.plt <- (yy.ma.plt + Vweight[ii] * y.expect.plt)
                  }
                  else cat("\n\nATTENTION:   NAs in model average response \n\n")
                }
                if (!cont) {
                  bg.resp <- unique(yy.ma[xx.interpol == 0])
                  if (ces.ans == 1) 
                    ces.abs <- 0.5
                  if (ces.ans == 2) 
                    ces.abs <- CES + bg.resp
                  if (ces.ans == 3) 
                    ces.abs <- CES * (1 - bg.resp) + bg.resp
                }
                if (cont) {
                  yy.ma <- exp(yy.ma)
                  yy.ma.plt <- exp(yy.ma.plt)
                  bg.resp <- unique(yy.ma[xx.interpol == 0])
                  ces.abs <- (CES + 1) * bg.resp
                }
                CED.lst <- approx(yy.ma, xx.interpol, xout = ces.abs)
                if (0) {
                  print("f.ced.ma")
                  print("y ma")
                  print(yy.ma[c(1, 100, 500)])
                  print("xx.inter")
                  print(xx.interpol[c(1, 100, 500)])
                  print("CED")
                  print(CED.lst$y)
                  print(bg.resp)
                }
                Vced.ma <- c(Vced.ma, CED.lst$y)
                if (!cont) {
                  if (dtype != 2) {
                    if (plot.type %in% c(2, 4, 6)) 
                      xline.plt <- log10(xline)
                    else xline.plt <- xline
                    if (plot.type %in% c(3, 4)) 
                      yy.ma.plt <- log10(yy.ma.plt)
                    if (plot.type == 5) 
                      yy.ma.plt <- asin(sqrt((yy.ma.plt)))
                    if (plot.type %in% c(1, 3, 5)) 
                      Vced.ma.plt <- Vced.ma
                    if (plot.type %in% c(2, 4, 6)) 
                      Vced.ma.plt <- log10(Vced.ma)
                    if (display_plots) {
                      lines(xline.plt, yy.ma.plt, lty = 2, col = group)
                      points(log10(Vced.ma[group]), -0.03, col = group)
                    }
                  }
                }
                if (cont) {
                  if (plot.type %in% c(1, 3)) 
                    xline.plt <- xline
                  if (plot.type %in% c(2, 4, 8)) 
                    xline.plt <- log10(xline)
                  if (plot.type %in% c(6, 8)) 
                    xline.plt <- sqrt(xline)
                  if (plot.type %in% c(3, 4, 6)) 
                    yy.ma.plt <- log10(yy.ma.plt)
                  if (plot.type %in% 7:8) 
                    yy.ma.plt <- sqrt(yy.ma.plt)
                  if (plot.type %in% c(1, 3)) 
                    Vced.ma.plt <- Vced.ma
                  if (plot.type %in% c(2, 4, 8)) 
                    Vced.ma.plt <- log10(Vced.ma)
                  if (plot.type %in% c(6, 8)) 
                    Vced.ma.plt <- sqrt(Vced.ma)
                  if (display_plots) {
                    lines(xline.plt, yy.ma.plt, lty = 2, col = group)
                    points(Vced.ma.plt[group], y.min, col = group)
                  }
                }
            }
            if (0) {
                x.tmp <- x
                x.tmp[x == 0] <- min(x[x > 0])/2
                x.min <- min(x.tmp)
                if (display_plots) {
                  f.graph.window(1)
                }
                plot(log10(xline), yy.ma.plt, xlim = log10(c(x.min, 
                  max(x))), type = "n")
                lines(log10(xline), yy.ma.plt)
            }
            if (0) {
                ans.all.tmp <- ans.all
                ans.all.tmp$plot.type <- 1
                ans.all.tmp$y <- ans.all$kk/ans.all$nn
                f.plot.frq(ans.all.tmp)
                lines((xx), yy.ma)
            }
            if (exists("track")) 
                print("f.ced.ma:   END 2")
            return(signif(Vced.ma, 4))
        }
    })
}



f.xx.ma <- function(ans.all, group, length.xx) {
    if (exists("track")) 
        print("f.xx.ma")
    ans.all$group <- group
    with(ans.all, {
        Vced <- numeric()
        Vced <- numeric()
        nr.gr <- max(covariate)
        for (ii in 1:nr.models) {
            if (!cont && dtype != 3) {
                CED.matr <- ans.all[[model.list[ii]]]$CED.matr
                if (length(CED.matr[, 1]) == 1 && nr.gr > 1) 
                  CED.matr <- matrix(CED.matr, nrow = nr.gr, 
                    ncol = 1, byrow = TRUE)
                Vced[ii] <- CED.matr[group, 1]
            }
            if (cont) {
                Vced[ii] <- ans.all[[model.list[ii]]]$CED[group]
            }
        }
        if (any(!is.na(Vced))) 
            xx.min <- min(Vced, na.rm = TRUE)/2
        else xx.min <- min(x[x > 0])/10
        Vced[Vced == Inf | Vced == -Inf] <- NA
        if (any(!is.na(Vced))) 
            xx.max <- max(Vced, na.rm = TRUE) * 2
        else xx.max <- max(x)
        if (xx.min == 0) 
            xx.min <- min(x[x > 0])/10
        if (xx.min == -Inf) 
            xx.min <- min(x[x > 0])/10
        xx <- exp(seq(from = logb(xx.min), to = logb(xx.max), 
            length = length.xx))
        if (min(xx) != 0) {
            xx <- c(0, xx)
        }
        if (exists("track")) 
            print("f.xx.ma:  END")
        return(xx)
    })
}



f.prepare.boot <- function(ans.all) {
    if (exists("track")) 
        print("f.prepare.boot")
    ans.all.bt <- ans.all
    with(ans.all, {
        if (dtype %in% c(10, 15, 250, 260)) {
            ans.all.bt$fct1 <- numeric()
            ans.all.bt$fct2 <- numeric()
            ans.all.bt$fct3 <- numeric()
            ans.all.bt$fct4 <- numeric()
            ans.all.bt$fct5 <- numeric()
            for (ii in 1:length(x)) {
                ans.all.bt$fct1 <- c(ans.all.bt$fct1, rep(fct1[ii], 
                  nn[ii]))
                ans.all.bt$fct2 <- c(ans.all.bt$fct2, rep(fct2[ii], 
                  nn[ii]))
                ans.all.bt$fct3 <- c(ans.all.bt$fct3, rep(fct3[ii], 
                  nn[ii]))
                ans.all.bt$fct4 <- c(ans.all.bt$fct4, rep(fct4[ii], 
                  nn[ii]))
                ans.all.bt$fct5 <- c(ans.all.bt$fct5, rep(fct5[ii], 
                  nn[ii]))
            }
            if (!do.MA) {
                ans.all.bt$x <- numeric()
                ans.all.bt$y.true <- numeric()
                for (ii in 1:length(x)) {
                  ans.all.bt$y.true <- c(ans.all.bt$y.true, rep(y.true[ii], 
                    nn[ii]))
                  ans.all.bt$x <- c(ans.all.bt$x, rep(x[ii], 
                    nn[ii]))
                }
            }
        }
        if (do.MA) {
            ans.all.bt$x <- numeric()
            ans.all.bt$y.true <- numeric()
            ans.all.bt$covariate <- numeric()
            for (ii in 1:length(xx.ma)) {
                ans.all.bt$y.true <- c(ans.all.bt$y.true, rep(yy.ma[ii], 
                  nn.ma[ii]))
                ans.all.bt$x <- c(ans.all.bt$x, rep(xx.ma[ii], 
                  nn.ma[ii]))
                ans.all.bt$covariate <- c(ans.all.bt$covariate, 
                  rep(covar.ma[ii], nn.ma[ii]))
            }
        }
        if (dtype %in% c(1, 10, 25, 26, 250, 260)) 
            sigma.bt <- rep(0, length(ans.all.bt$x))
        if (!do.MA) {
            for (jj in (1:max(fct3))) if (dtype %in% c(1, 10, 
                25, 26, 250, 260)) 
                sigma.bt <- sigma.bt + sqrt(MLE[jj]) * (ans.all.bt$fct3 == 
                  jj)
            if (dtype %in% c(5, 15)) {
                inter.sigma.bt <- rep(sqrt(inter.var), length(ans.all.bt$x))
                intra.sigma.bt <- rep(sqrt(intra.var), length(ans.all.bt$x))
            }
        }
        if (do.MA) {
            for (jj in (1:max(covariate))) {
                if (dtype %in% c(1, 10, 25, 26, 250, 260)) 
                  sigma.bt <- sigma.bt + sigma.ma[jj] * (ans.all.bt$covariate == 
                    jj)
                if (dtype %in% c(5, 15)) {
                  inter.sigma.bt <- rep(inter.sigma.ma, length(ans.all.bt$x))
                  intra.sigma.bt <- rep(intra.sigma.ma, length(ans.all.bt$x))
                }
            }
        }
        if (dtype %in% c(1, 10, 25, 26, 250, 260)) 
            ans.all.bt$sigma.bt <- sigma.bt
        if (dtype %in% c(5, 15)) {
            ans.all.bt$inter.sigma.bt <- inter.sigma.bt
            ans.all.bt$intra.sigma.bt <- intra.sigma.bt
        }
        if (!do.MA) 
            ans.all.bt$par.start <- MLE
        if (dtype == 10) 
            ans.all.bt$dtype <- 1
        if (dtype == 15) 
            ans.all.bt$dtype <- 5
        if (exists("track")) 
            print("f.prepare.boot  END")
        return(ans.all.bt)
    })
}



f.gener.con <- function(ans.all) {
    if (exists("track")) 
        print("f.gener.con")
    with(ans.all, {
        noise <- numeric()
        df.corr <- length(x) - nrp
        if (dtype %in% c(5, 15)) {
            inter.sigma.bt <- inter.sigma.bt * sqrt(df1/rchisq(1, 
                df1))
            intra.sigma.bt <- intra.sigma.bt * sqrt(df2/rchisq(1, 
                df2))
            eps <- rnorm(no.of.nests, 0, inter.sigma.bt)
            for (ii in 1:no.of.nests) {
                noise <- c(noise, rnorm(nest.size[ii], eps[ii], 
                  intra.sigma.bt))
                ans.all$dtype <- 1
            }
        }
        else {
            sigma.bt <- sigma.bt * sqrt(df.corr/rchisq(1, df.corr))
            noise <- rnorm(length(x), 0, sigma.bt)
        }
        if (dtype == 25 | dtype == 250) 
            y <- y.true + (noise)
        else if (dtype == 26 | dtype == 260) 
            y <- (sqrt(y.true) + (noise))^2
        else y <- y.true * exp(noise)
        if (exists("track")) 
            print("f.gener.con:  END")
        return(y)
    })
}



f.show.ma <- function(ans.all) {
    if (exists("track")) 
        print("f.show.ma")
    with(ans.all, {
        nl <- "\n"
        sp <- " "
        show <- paste("\nversion:", ans.all$PRversion)
        show <- paste(show, nl, "model averaging results")
        show <- paste(show, nl, "dtype", sp, dtype)
        show <- paste(show, nl, "selected", sp, ans.all$select.name)
        if (ans.all$remove.name != "none") 
            show <- paste(show, nl, "removed", sp, ans.all$remove.name)
        if (cond.ans != 3) 
            show <- paste(show, "\n fit cond:", cond.ans)
        show <- paste(show, factor.name)
        show <- paste(show, nl, "dose scaling:", sf.x)
        show <- paste(show, nl, "conf level:", conf.lev)
        show <- paste(show, nl, "number of runs:", nr.boot.ma)
        if (!cont) {
            if (CED.model) 
                switch(ces.ans, show <- paste(show, nl, "CED = ED50"), 
                  show <- paste(show, nl, "addit. risk", sp, 
                    CES), show <- paste(show, nl, "extra risk", 
                    sp, CES))
            nr.gr <- length(MA$conf.int[, 1])
        }
        if (cont) 
            show <- paste(show, nl, "CES", sp, CES)
        if (!is.na(MA$Vlower.ma[1])) {
            nr.gr <- length(MA$Vlower.ma)
            show <- paste(show, nl, "BMD CI")
            for (ii in 1:nr.gr) {
                show <- paste(show, nl, signif(MA$conf.int[ii, 
                  2], 2), signif(MA$conf.int[ii, 3]), sep = "  ")
            }
        }
        if (exists("track")) 
            print("f.show.ma: END")
        return(show)
    })
}



f.mtext <- function(ans.all, display_plots = TRUE) {
    if (exists("track")) 
        print("f.mtext")
    modelname <- ans.all$modelname
    show <- ans.all$show
    if (ans.all$show != "") {
        if (length(ans.all$y.lim.plt) > 0) 
            y.lim <- ans.all$y.lim.plt
        else y.lim <- ans.all$xy.lim[4:5]
        loc.txt <- y.lim[2] + 0.05 * (y.lim[2] - y.lim[1])
        if (display_plots == TRUE) {
            mtext(show, 4, 1, adj = 0, padj = 1, las = 1, at = loc.txt, 
                cex = 0.7)
        }
        if (ans.all$model.ans %in% c(31, 32, 33, 45, 55)) 
            cexmain <- 0.8
        else cexmain <- 1.2
        if (ans.all$output)
            if (display_plots == TRUE) {
                title(main = paste("\n\n", modelname), cex.main = cexmain, 
                    font.main = 1)
            }
    }
    else if (!ans.all$fitted && ans.all$heading == "") {
        title(main = "\nstarting values")
    }
    if (exists("track")) 
        print("f.mtext:  END")
    return(invisible())
}



f.boot.ma <- function(ans.all, interactive_mode = TRUE, results_env = NULL, display_plots = TRUE, filename = NULL, output_type = NULL) {
    if (exists("track")) 
        print("f.boot.ma")
    if (ans.all$seed.bt != 0) 
        set.seed(ans.all$seed.bt)
    date.0 <- date()
    with(ans.all, {
        nr.gr <- max(ans.all$covariate)
        if (dtype != 2) {
            dum.contr <- xy.lim[1]/10
            xline <- exp(seq(from = logb(2 * dum.contr), to = logb(xy.lim[3]), 
                length = 1000))
            name.wapp <- paste("c", yans, "ma.plot", sep = "")
            if (display_plots) {
                f.graph.window(1, WAPP = WAPP, title = "bootstrap curves based on model averaging", 
                name.wapp = name.wapp, plotprefix = ans.all$plotprefix, 
                svg.plots = svg.plots, output = output, filename = filename,
            output_type = output_type)
            }
            if (!cont) 
                ans.all$xy.lim[4:5] <- c(0, 1)
            ans.all$xy.lim[1] <- dum.contr
            ans.all$heading <- ""
            if (covar.no == 0) 
                ans.all$color <- ans.all$color[-1]
            if (!cont) 
                f.plot.frq(ans.all)
            if (cont) {
                ans.all$fct2 <- covariate
                ans.all$plt.mns <- 3
                ans.all$heading <- ""
                ans.all$y.min <- f.plot.con(ans.all, display_plots = display_plots)$y.lim.plt[1]
            }
            if (display_plots) {
                title(main = "bootstrap curves \nbased on model averaging", 
                cex.main = 1)
            }
        }
        else xline <- NA
        ans.all$MA <- f.ced.ma(ans.all, xline = xline, first.call = TRUE, display_plots = display_plots)
        cat("\n The weights used in model averaging are:\n")
        print(ans.all$MA$Vweight)
        date.0 <- date()
        cat("\nStart of MA bootstrap runs ...\n")
        ans.all.bt <- ans.all
        ans.all.bt$output <- FALSE
        if (dtype %in% c(2, 4, 6)) {
            yy.true <- ans.all$MA$yy.ma
            ans.all.bt$x <- ans.all$MA$xx.ma
            ans.all.bt$nn <- ans.all$MA$nn.ma
            ans.all.bt$covariate <- ans.all$MA$covar.ma
        }
        if (cont) {
            ans.all.bt$xx.ma <- ans.all$MA$xx.ma
            ans.all.bt$yy.ma <- ans.all$MA$yy.ma
            ans.all.bt$nn.ma <- as.numeric(ans.all$MA$nn.ma)
            ans.all.bt$covar.ma <- ans.all$MA$covar.ma
            if (dtype %in% c(1, 10, 25, 26, 250, 260)) 
                ans.all.bt$sigma.ma <- ans.all$MA$sigma.ma
            if (dtype %in% c(5, 15)) {
                ans.all.bt$inter.sigma.ma <- ans.all$MA$inter.sigma.ma
                ans.all.bt$intra.sigma.ma <- ans.all$MA$intra.sigma.ma
            }
            ans.all.bt <- f.prepare.boot(ans.all.bt)
        }
        ans.all.bt$gui <- TRUE
        ans.all.bt$no.CI <- TRUE
        ans.all.bt$full.ans <- 2
        Vced.ma <- numeric()
        ced.ma.matr <- matrix(nrow = nr.boot.ma, ncol = nr.gr)
        if (!WAPP && output) 
            cat("run ")
        for (ii in (1:nr.boot.ma)) {
            if (!WAPP && output) 
                cat(" ", ii)
            else cat("run", ii, "\n")
            if (efsa.tool) 
                shiny::incProgress(1/nr.boot.ma, detail = paste("bootstrap", 
                  ii, "out of", nr.boot.ma))
            ans.all.bt$dtype <- ans.all$dtype
            if (!cont) {
                yy.true.tmp <- yy.true
                if (dtype == 6) {
                  beta <- (alfa.mle * (1 - yy.true))/yy.true
                  yy.true.tmp <- rbeta(length(x), alfa.mle, beta)
                }
                ans.all.bt$kk <- rbinom(length(x), ans.all$MA$nn, 
                  yy.true.tmp)
                ans.all.bt$y <- ans.all.bt$kk/ans.all$MA$nn
            }
            if (cont) {
                ans.all.bt$y <- f.gener.con(ans.all.bt)
                ans.all.bt$yy <- ans.all.bt$y
            }
            if (dtype == 2) 
                ans.all.bt$y <- ans.all.bt$kk
            if (cont) 
                if (any(c(!ans.all$EXP$trend, !ans.all$HILL$trend, 
                  !ans.all$INVEXP$trend, !ans.all$LOG$trend))) 
                  ans.all.bt$quick.ans <- 4
            if (!cont) 
                ans.all.bt <- f.select.cat(ans.all.bt, output = FALSE, interactive_mode = interactive_mode)
            if (cont) {
                ans.all.bt$fct4 <- rep(1, length(ans.all.bt$x))
                ans.all.bt$dtype <- 1
                ans.all.bt$dtype.0 <- 1
                ans.all.bt$output <- FALSE
                ans.all.bt <- f.select.con(ans.all.bt, interactive_mode = interactive_mode, results_env = results_env, display_plots = display_plots)
            }
            Vced.ma <- f.ced.ma(ans.all.bt, xline = xline, first.call = FALSE, display_plots = display_plots)
            ced.ma.matr[ii, ] <- Vced.ma
        }
        Vlower.ma <- numeric()
        Vupper.ma <- numeric()
        for (jj in 1:nr.gr) {
            Vlower.ma[jj] <- signif(quantile(ced.ma.matr[, jj], 
                (1 - conf.lev)/2, na.rm = TRUE), 3)
            Vupper.ma[jj] <- signif(quantile(ced.ma.matr[, jj], 
                (1 + conf.lev)/2, na.rm = TRUE), 3)
        }
        ans.all$MA$nr.boot.ma <- nr.boot.ma
        ans.all$MA$Vlower.ma <- Vlower.ma
        ans.all$MA$Vupper.ma <- Vupper.ma
        ans.all$MA$conf.lev <- conf.lev
        ans.all$MA$conf.int.ma <- data.frame(subgroup = covar.txt, 
            BMDlower.ma = Vlower.ma, BMDupper.ma = Vupper.ma)
        if (nr.gr == 1) 
            ans.all$MA$conf.int.ma <- data.frame(subgroup = "all", 
                BMDlower.ma = Vlower.ma, BMDupper.ma = Vupper.ma)
        ans.all$MA$ced.ma.matr <- ced.ma.matr
        CI.row.ma <- numeric()
        for (ii in 1:nr.gr) CI.row.ma <- c(CI.row.ma, Vlower.ma[ii], 
            Vupper.ma[ii])
        ans.all$MA$CI.row.ma <- CI.row.ma
        if (dtype != 2) {
            ans.all$show <- f.show.ma(ans.all)
            ans.all$output <- FALSE
            f.mtext(ans.all, display_plots = display_plots)
            ans.all$output <- TRUE
        }
        if (!ans.all$gui) {
            cat("\n\nduration of bootstrap calculations:\n")
            print(date.0)
            print(date())
        }
        if (any(is.na(ced.ma.matr))) 
            cat("\n\nATTENTION: \n        There are NAs in the vector of bootstrap CEDs, this indicates a problem in interpolation\n\n")
        if (WAPP) 
            dev.off()
        if (exists("track")) 
            print("f.boot.ma: END")
        return(ans.all)
    })
}



f.store.results <- function(ans.all, store.name = 0, add.ext = 0, interactive_mode = NULL) {
    if (exists("track")) 
        print("f.store.results")
    if (ans.all$quick.ans == 1) {
        ans.all$EXP <- NULL
        ans.all$HILL <- NULL
    }
    if (interactive_mode == FALSE) {
        return(ans.all)
    } else {
        if (store.name == 0) 
            ans.all$store.name <- eval(as.character(parse(prompt = "\n\ngive name for file to store results (or type 0 if none)  >   ")))
            ans.all$date <- date()
            if (ans.all$store.name != 0) {
                f.assign(eval(ans.all$store.name), ans.all)
            }
            return(ans.all$store.name)
    }
}



f.version <- function() {
    PRversion <- "71.1"
    PRversion
}



f.control <- function(lev = 3) {
    lst <- list()
    switch(lev, {
        lst$eval.max <- 50
        lst$iter.max <- 40
        lst$rel.tol <- 0.001
        lst$x.tol <- 0.015
        lst$step.min <- 0.00022
    }, {
        lst$eval.max <- 100
        lst$iter.max <- 75
        lst$rel.tol <- 1e-06
        lst$x.tol <- 0.00015
        lst$step.min <- 2.2e-07
    }, {
        lst$eval.max <- 1000
        lst$iter.max <- 750
    })
    return(lst)
}



f.con <- function(ans.all, list.logic = F, indep_var_choice = NULL, Vyans_input = NULL, interactive_mode = TRUE, covariates = NULL, custom_CES = NULL, model_selection = NULL, lower_dd = NULL, upper_dd = NULL, selected_model = NULL, adjust_CES_to_group_SD = NULL, model_averaging = NULL, num_bootstraps = NULL, results_env = NULL, display_plots = TRUE) {
    if (exists("track")) 
        print("f.con")
    f.assign(".Pr.last", ans.all)
    if (list.logic) {
        cat("\n You have chosen previous results concerning:   ")
        cat(ans.all$odt.name)
        PRversion <- ans.all$PRversion
        if (ans.all$model.fam > 0) {
            cat("\nfitted model: ")
            cat(ans.all$modelname, "\n")
            if (length(ans.all$loglik) == 1) 
                cat("loglik:", ans.all$loglik, "\n")
            cat(ans.all$show)
        }
        if (ans.all$boot) {
            cat("\nBootstrap data are available")
            nruns <- length(ans.all$CED.boot[, 1])
            cat("\nnruns:", nruns, "\n")
        }
        if (length(ans.all$SINGMOD) > 0) 
            ans.all <- f.move.sublist(ans.all, ans.all$SINGMOD)
        if (ans.all$model.ans %in% 41:42 && ans.all$nr.cc > 1) {
            with(ans.all, {
                npar <- length(MLE)
                Vcc <- MLE[(npar - nr.dd + 1):npar]
                Vsd <- sqrt(MLE[1:nr.var])
                if (display_plots) {
                  f.graph.window(1)
                  plot(log10(Vcc), log10(Vsd), xlab = "log.sd", 
                  ylab = "log.c")
                }
                f.press.key.to.continue()
            })
        }
        if (length(ans.all$xans) > 1) 
            ans.all$ans.plt <- 0
        if (ans.all$nr.cc > 1 && ans.all$incr.decr.no == 0) {
            ans.all <- with(ans.all, {
                npar <- length(MLE)
                Vcc <- MLE[(nr.var + nr.aa + nr.bb + 1):(nr.var + 
                  nr.aa + nr.bb + nr.cc)]
                Vsd <- sqrt(MLE[1:nr.var])
                if (length(Vsd) == length(Vcc)) {
                  gr.txt <- f.pars(ans.all)$gr.txt
                  gr.lab <- letters[1:length(gr.txt)]
                  show <- "subgroups"
                  for (i in 1:length(gr.txt)) show <- paste(show, 
                    nl, gr.lab[i], tb, gr.txt[i])
                  if (display_plots) {
                    f.graph.window(1)
                  }
                  plot(Vsd, log(Vcc), pch = as.character(gr.lab))
                  abline(0, 7)
                  loc.txt <- max(log(Vcc))
                  mtext(show, 4, 1, adj = 0, padj = 1, las = 1, 
                    at = loc.txt, cex = 0.7)
                  cat("\n\n\nSee plot for relationship between c and sd\n\n")
                  f.press.key.to.continue()
                  if (model.ans == 6 && fitted) 
                    ans.all <- f.fit.model6(ans.all)
                }
                return(ans.all)
            })
        }
        ans.all$par.start <- ans.all$MLE
        if (ans.all$plot.type > 8) 
            ans.all$plot.type <- 2
        ans.all$PRversion <- ans.all$version.old
        if (ans.all$model.fam > 0) 
            f.plot.all(ans.all)
        else f.plot.gui(ans.all)
        ans.all$PRversion <- PRversion
        if (ans.all$dtype != 3 && length(ans.all$gr.txt) > 1) 
            f.explain.marks(ans.all)
        if (length(ans.all$EXP$conf.int[, 1]) > 1) {
            f.press.key.to.continue()
            cat("\n\nplot of BMD CIs for each subgroup\n\n")
            f.plot.CI(ans.all)
        }
    }
    main.ans <- 2
    if (list.logic) 
        main.ans <- 12
    if (ans.all$model.ans > 0) {
        if (!ans.all$fitted) 
            ans.all <- f.start.con(ans.all)
        main.ans <- 14
    }
    while (T) {
        switch(main.ans, {
            ans.all$change <- rep(F, ans.all$nrQ)
            ans.all <- f.change.settings(ans.all, choose = T)
            if (ans.all$fitted) ans.all$show <- f.show.con(ans.all)
            if (!ans.all$fitted) {
                if (ans.all$model.ans %in% 38:40) ans.all <- f.quick.con(ans.all, indep_var_choice = indep_var_choice, Vyans_input = Vyans_input, covariates = covariates, custom_CES = custom_CES, model_selection = model_selection, lower_dd = lower_dd, upper_dd = upper_dd, interactive_mode = interactive_mode, adjust_CES_to_group_SD = adjust_CES_to_group_SD, model_averaging = model_averaging, num_bootstraps = num_bootstraps, results_env = results_env, display_plots = display_plots) else {
                  ans.all <- f.execute(ans.all)
                  if (ans.all$ans.m6.sd != 2) ans.all <- f.start.con(ans.all, 
                    adjust = F, fitted = F, tmp.quick = F)
                }
            } else {
                ans.all <- f.execute(ans.all)
                f.plot.all(ans.all)
            }
            model.names <- f.expect.con(name = TRUE)
            ans.all$modelname <- model.names[ans.all$model.ans]
        }, {
            ans.all <- f.clear(ans.all)
            if (ans.all$quick.ans == 1) {
                ans.all <- f.choose.model(ans.all)
                if (ans.all$model.ans %in% c(48, 50)) {
                  ans.all$change[24] <- T
                  ans.all <- f.change.settings(ans.all)
                  ans.all <- f.execute(ans.all)
                }
                if (ans.all$model.ans == 6) ans.all$ans.m6.sd <- 1
                ans.all$new.model <- T
            }
            if (ans.all$quick.ans == 1) {
                ans.all$modelname <- ans.all$model.names[ans.all$model.ans]
                ans.all <- f.start.con(ans.all)
                if (ans.all$adjust.start) {
                  ans.all <- f.start.con(ans.all, adjust = T, 
                    fitted = T, tmp.quick = F)
                }
            }
            if (ans.all$quick.ans > 1) {
                ans.all <- f.quick.con(ans.all, indep_var_choice = indep_var_choice, Vyans_input = Vyans_input, covariates = covariates, custom_CES = custom_CES, model_selection = model_selection, lower_dd = lower_dd, upper_dd = upper_dd, interactive_mode = interactive_mode, adjust_CES_to_group_SD = adjust_CES_to_group_SD, model_averaging = model_averaging, num_bootstraps = num_bootstraps, results_env = results_env, display_plots = display_plots)
                
                # Human-readable list of models
                list.of.models <- c("exponential", "Hill", "inverse exponential", "lognormal DR")
      
                if (interactive_mode == FALSE && !is.null(selected_model)) {
                    # Find the index of the selected model in the list of models and proceed
                  # TODO Error checking on the user input (selected_model parameter must be one of the list.of.models)
                    model.fam <- which(tolower(list.of.models) == tolower(selected_model)) # assuming selected_model is the name of the model
                    ans.all$model.fam <- model.fam
                } else {
                    # In case of interactive mode or no selected model, show a menu to continue
                    model.fam <- menu(list.of.models, title = "\nWhich model do you want to continue with?  ")
                    ans.all$model.fam <- model.fam
                }
                
                # # Use the model.fam to determine further steps
                # ans.all <- switch(ans.all$model.fam,
                #     # Assuming EXP, HILL, INVEXP, LOGN match indices of list.of.models
                #     f.move.sublist(ans.all, ans.all$EXP),
                #     f.move.sublist(ans.all, ans.all$HILL),
                #     f.move.sublist(ans.all, ans.all$INVEXP),
                #     f.move.sublist(ans.all, ans.all$LOGN),
                #     stop("Invalid model family selection.")
                # )
                # # Store results and other tasks
                # f.store.results(ans.all)
                
                # # If 'selected_model' was used, 'model.fam' will already be set, 
                # # otherwise, it will be set by the menu choice
                # ans.all$modelname <- list.of.models[ans.all$model.fam]
                
                # # Quick analysis reset to step into the next interaction until an "exit condition" is met
                # ans.all$quick.ans <- 1
            
                # # Set ans.plt to 0 to probably reset the plot or interaction
                # ans.all$ans.plt <- 0
              
              
              
              
              
              # list.of.models <- c("exponential", "Hill", "inverse exponential", 
              #     "lognormal DR")
              #   list.of.models <- list.of.models[1:ans.all$nr.models]
              #   if (length(ans.all$Vyans) >= 1 && length(ans.all$HILL) > 
              #     0) model.fam <- menu(list.of.models, title = "\nWhich model do you want to continue with?  ") else model.fam <- 1
                if (ans.all$quick.ans < 6) switch(model.fam, 
                  ans.all <- f.move.sublist(ans.all, ans.all$EXP), 
                  ans.all <- f.move.sublist(ans.all, ans.all$HILL), 
                  ans.all <- f.move.sublist(ans.all, ans.all$INVEXP), 
                  ans.all <- f.move.sublist(ans.all, ans.all$LOGN)) else ans.all <- f.move.sublist(ans.all, 
                  ans.all$SINGMOD)
                f.store.results(ans.all, interactive_mode = interactive_mode)
                if (0) switch(model.fam, ans.all$heading <- "Exponential model", 
                  ans.all$heading <- "Hill model", ans.all$heading <- "Inverse exponential model", 
                  ans.all$heading <- "Lognormal DR model")
                ans.all$model.fam <- model.fam
                ans.all$quick.ans <- 1
                ans.all$modelname <- ans.all$model.names[ans.all$model.ans]
                ans.all$ans.plt <- 0




            }
          

          
        }, {
            if (ans.all$fitted & !ans.all$new.model) ans.all$par.start <- ans.all$MLE
            if (0) {
                ans.dd <- 2
                if (max(ans.all$fct5) > 1) ans.dd <- menu(c("yes", 
                  "no"), title = "\nDo you want to use previous MLE from model with one d as start value?\n")
                if (ans.dd == 1) {
                  nr.dd <- max(ans.all$fct5)
                  dd.start <- ans.all$MLE[length(ans.all$MLE)]
                  ans.all$par.start <- c(ans.all$MLE, rep(dd.start, 
                    nr.dd - 1))
                  ans.all$scale.dum <- c(ans.all$scale.dum, rep(1, 
                    nr.dd - 1))
                }
            }
            ans.all <- f.start.con(ans.all, adjust = T, fitted = ans.all$fitted, 
                tmp.quick = F)
        }, {
            ans.all <- f.mm4.con(ans.all)
            while (ans.all$out.ans == 1) ans.all <- f.mm4.con(ans.all)
        }, {
            ans.all <- f.mm5.con(ans.all)
        }, {
            if (ans.all$fitted == F) cat("\nFirst fit the model using option 4\n") else if (ans.all$model.ans == 
                11) cat("\nCED not defined for full model\n") else if (ans.all$model.ans == 
                1) cat("\nCED not defined for null model\n") else ans.all <- f.mm6.con(ans.all)
        }, {
            if (is.na(ans.all$CED[1]) && ans.all$model.ans != 
                59) cat("\nYou did not calculate CED! First use option 6 from main menu\n") else {
                ans.all$nruns <- eval(parse(prompt = "\ngive the number of bootstrap runs you want to be done > "))
                ans.all$plot.ans <- 1
                ans.all$plot.ans <- menu(c("no", "yes (reduces speed)"), 
                  title = "\ndo you want plots of each sampling run?")
                ans.all <- f.mm7.con(ans.all)
            }
        }, {
            ans.all <- f.mm8.con(ans.all)
        }, {
            ans.all <- f.mm10.con(ans.all)
        }, {
            with(ans.all, {
                if (model.ans != 2) {
                  cat("\nonly model 2 implemented\n")
                  return(invisible())
                }
                regr.par.matr <- f.pars(ans.all)$regr.par.matr
                for (ii in 1:length(regr.par.matr[, 1])) {
                  aa <- logb(regr.par.matr[ii, 1])
                  bb <- regr.par.matr[ii, 2]
                  if (max(fct3) == 1) SD <- sqrt(MLE[1]) else SD <- sqrt(MLE[ii])
                  ADI <- eval(parse(prompt = paste("give value for ADI (or MRL)", 
                    "  > ")))
                  kkk <- eval(parse(prompt = paste("give value for k, in tol.lim = x + k*sd", 
                    "  > ")))
                  WP <- (logb(ADI) - kkk * SD - aa)/bb
                  cat("\n\nThe withdrawal period is:", WP, "\n")
                }
            })
        }, {
            with(ans.all, {
                cat(paste(1:npar, ": ", text.par[1:npar], "\n"))
                cat("Give number(s) of the (consecutive) parameter(s)\n")
                ans.all$group <- eval(parse(prompt = paste(" -------- > ")))
                MLE <- ans.all$MLE
                if (length(MLE) > 100 && !WAPP) {
                  cat("\n\nnumber of parameters exceeds 100\n")
                  cat("Consider to relax the convergence criteria\n\n")
                  ans.all$change[23] <- T
                  ans.all$change[nrQ] <- F
                  ans.all <- f.change.settings(ans.all)
                }
                conf.int.0 <- ans.all$conf.int
                ans.all$trace <- T
                ans.all <- f.CI(ans.all)
                name.ci <- text.par[ans.all$group]
                confint.spec <- ans.all$conf.int
                if (ans.all$sf.x != 1) {
                  cat("\nATTENTION: the scaling factor on dose of", 
                    ans.all$sf.x, "was not taken into account here")
                  cat("\nwhich should be noted when one of confidence intervals relate to a CED\n")
                }
                dimnames(confint.spec) <- list(name.ci, NULL)
                ans.all$conf.int <- conf.int.0
                ans.all$confint.spec <- confint.spec
                cat("\n\n the confidence intervals are:\n")
                print(confint.spec)
                f.store.results(ans.all, interactive_mode = interactive_mode)
            })
        }, {
            return(T)
        })
      if (interactive_mode == FALSE) {
        main.ans <- 12
      } else {
        cat("\n\n--------\nMAIN MENU :\n--------\n\n")
        main.ans <- menu(c("Change settings", "Choose (another) model", 
            "Choose other startvalues", "Fit model", "Plot results", 
            "Calculate CED: point estimate or confidence interval", 
            "Generate bootstrap runs based on fitted model", 
            "Calculate CED distribution for the animal", "Calculate effect size for given exposure in the animal", 
            "Calculate withdrawal period for veterinary medicines", 
            "Calculate confidence interval for model parameter(s)", 
            "End session"), title = "What do you want to do ? \n")
      }
      }
}



f.nr.replicates <- function(ans.all) {
    if (exists("track")) 
        print("f.nr.replicates")
    with(ans.all, {
        if (length(xans) > 1) {
            x <- data.0[, xans[1]]
        }
        if (quick.ans == 1) 
            if (length(fct3) > 1) 
                dum.nn <- table(x, fct1, fct2, fct3, fct5)
            else dum.nn <- table(x, fct1, fct2, fct5)
        if (quick.ans > 1) 
            dum.nn <- table(x, covariate)
        if (dtype %in% c(4, 6, 10, 15, 250, 260)) 
            dum.nn <- nn
        dum.nn <- as.numeric(dum.nn)
        dum.nn <- dum.nn[dum.nn != 0]
        return(dum.nn)
    })
}



f.quick.con <- function(ans.all,
                        indep_var_choice = NULL,
                        Vyans_input = NULL,
                        covariates = NULL,
                        custom_CES = NULL,
                        model_selection = NULL,
                        interactive_mode = TRUE,
                        lower_dd = NULL,
                        upper_dd = NULL,
                        adjust_CES_to_group_SD = NULL,
                        model_averaging = NULL,
                        num_bootstraps = NULL,
                        results_env = NULL,
                        display_plots = TRUE) {
    if (exists("track")) 
        print("f.quick.con")
        message(paste0("indep_var_choice: ", indep_var_choice))
    if (ans.all$WAPP) {
        ans.all$gui <- TRUE
        ans.all$interrupt <- FALSE
    }
    ans.all <- with(ans.all, {
        output <- ans.all$output
        if (gui) {
            if (ans.all$quick.ans != 6) 
                ans.all$quick.ans <- 3
            if (output) 
                if (!WAPP && !simu) {
                  cat("\nclick in this window to see current output\n")
                  f.press.key.to.continue()
                }
            ans.all$data.0 <- odt
            f.overlap()
            ans.all$cont <- f.cont(ans.all)
            Vyans <- yans
            Vnans <- 0
            Vsans <- 0
            if (dtype %in% c(10, 15)) {
                Vsans <- sans
                Vnans <- nans
            }
            if (length(Vnans) == 1 && length(Vyans) > 1) 
                Vnans <- rep(Vnans, length(Vyans))
            ans.all$const.var <- FALSE
            version.old <- ans.all$PRversion
            if (version.old[1] != PRversion) 
                ans.all <- f.adjust.saved(ans.all)
            if (ans.all$WAPP && ans.all$nruns == 0 && dtype %in% 
                c(5, 15)) 
                ans.all$nruns <- 200
            if (ans.all$WAPP && ans.all$do.MA) 
                ans.all$nr.models <- 4
            if (length(Vyans) > 1) 
                ans.all$auto.detlim <- TRUE
            if (DA.ans == 2) 
                ans.all <- f.add.dosecols(ans.all)
        }
        if (!gui) {
            varnames <- ans.all$varnames
            nvar <- ans.all$nvar
            dtype <- ans.all$dtype
            odt <- ans.all$odt
            if (ans.all$xans[1] == 0) {
                ans.all$change[1] <- T
                ans.all <- f.change.settings(ans.all, indep_var_choice = indep_var_choice)
            }
            cat("\n")
            cat(paste(1:nvar, varnames[1:nvar], "\n"))
            cat("Give number(s) of the response(s) you want to analyse\n")
            if (!is.null(Vyans_input)) {
              if (is.numeric(Vyans_input) || all(Vyans_input %in% varnames)) {
                Vyans <- if (is.numeric(Vyans_input)) Vyans_input else match(Vyans_input, varnames)
                ans.all$Vyans <- Vyans
                zero.ind <- if (length(Vyans) > 1) any(data.0[, Vyans] == 0, na.rm = TRUE) else F
              } else {
                stop("Vyans_input must be a vector of either numeric column indices or names that exist in 'varnames'.")
              }
            } else {
            Vyans <- eval(parse(prompt = paste(" -------- > ")))
            }
            zero.ind <- F
            if (length(Vyans) > 1) 
                zero.ind <- any(data.0[, Vyans] == 0, na.rm = T)
            Vsans <- 0
            Vnans <- 0
            if (dtype %in% c(10, 15, 250)) {
                cat("\n")
                cat(paste(1:nvar, varnames[1:nvar], "\n"))
                cat("Give number(s) of the variation statistic associated to these endpoints\n")
                Vsans <- eval(parse(prompt = paste(" -------- > ")))
                if (length(Vsans) == 1 && length(Vyans) > 1) 
                  Vsans <- rep(Vsans, length(Vyans))
                sd.se <- menu(c("standard deviations", "standard errors"), 
                  title = "Do you have standard deviations or standard errors associated \n                          with the means?")
                cat("\n")
                cat(paste(1:nvar, varnames[1:nvar], "\n"))
                cat("Give number(s) of the group size(s) associated to these endpoints\n")
                Vnans <- eval(parse(prompt = paste(" -------- > ")))
                if (length(Vnans) == 1 && length(Vyans) > 1) 
                  Vnans <- rep(Vnans, length(Vyans))
                ans.all$Vsans <- Vsans
                ans.all$sd.se <- sd.se
                ans.all$Vnans <- Vnans
            }
            ans.all$covar.dd <- 0
            if (dtype == 5) 
                ans.all$nest.no <- menu(varnames[1:nvar], title = "Give number of the nested factor \n (type 0 to skip nested analysis)")
            if (ans.all$DA.ans == 1) {
                if (!is.null(covariates)) {
                  if (is.numeric(covariates) || covariates %in% varnames) {
                    covar.no <- if (is.numeric(covariates)) covariates else match(covariates, varnames)
                    ans.all$covar.no <- covar.no
                  } else {
                    stop("covariates must be either a numeric column index or a name that exists in 'varnames'.")
                  }
                } else if (interactive_mode && ans.all$DA.ans == 1) {
                  covar.no <- menu(varnames[1:nvar], title = "Give number of factor serving as potential covariate (e.g.sex)\n                     type 0 if none --- \n  ")
                  ans.all$covar.no <- covar.no
                }
            }
            if (0) 
                if (ans.all$quick.ans == 4 && covar.no > 0) {
                  covar.dd <- menu(c("no", "yes"), title = "Do you want parameter d to be covariate dependent? \n  ")
                  covar.dd <- covar.dd - 1
                  ans.all$covar.dd <- covar.dd
                  ans.all$CES <- -1
                  fix.ans <- 2
                  cc.fix <- NA
                  if (fix.ans == 1) 
                    ans.all$cc.fix <- eval(parse(prompt = "\n\ngive fixed value for parameter c >   "))
                }
            if (!ans.all$covar.dd) {
              if (interactive_mode == FALSE) {
                # Setting CES based on the parameter provided
                ans.all$NES.ans <- adjust_CES_to_group_SD
              } else {
                ans.all$NES.ans <- menu(c("no", "yes"), title = "Do you want to adjust CES to within group SD?\n")
              }
              if (ans.all$NES.ans == 1) {
                if (interactive_mode == FALSE ) {
                  ans.all$CES <- custom_CES
                } else {
                  ans.all$CES <- eval(parse(prompt = paste("\nGive value for CES (always positive)\n                           type 0 to avoid calculation of CIs  > ")))
                }
                ans.all$CES <- abs(ans.all$CES)
              }
            }
            if (quick.ans == 6) {
                fct2.fact <- data.0[, covar.no]
                levels.all <- levels(factor(fct2.fact))
                nr.lev <- length(levels.all)
                cat("\n: Give number associated with the reference\n")
                cat(" ")
                cat(paste(1:nr.lev, ":", levels.all, "\n"))
                ans.all$ref.lev <- eval(parse(prompt = paste(" -------- > ")))
            }
            if (length(Vyans) == 1) 
                interrupt <- TRUE
            if (length(Vyans) > 1) {
              if (interactive_mode == TRUE) {
                int.ans <- menu(c("no", "yes (opportunity to save results per endpoint)"), 
                title = "Do you want to interrupt calculations after each endpoint? \n")
              } else {
                int.ans <- 1
              }  
              if (int.ans == 2) 
                interrupt <- TRUE
              else interrupt <- FALSE
            }
            ans.all$interrupt <- interrupt
            ans.all$auto.detlim <- FALSE
            if (length(Vyans) > 1) {
                if (zero.ind && int.ans == 1) {
                  auto.ans <- menu(c("yes", "apply automatic detection limit"), 
                    title = "Do you want to provide detection limit for each separate response, or automically use smallest nonzero observation?\n")
                  switch(auto.ans, ans.all$auto.detlim <- FALSE, 
                    ans.all$auto.detlim <- TRUE)
                }
            }
            if ((ans.all$CES != 0 || ans.all$NES.ans == 2)) {
                if (ans.all$quick.ans == 6) 
                  ans.all$nruns <- eval(parse(prompt = "\nprovide number of bootstrap runs, or type 0 to get profile likelihood CI > "))
                else if (length(ans.all$xans) == 1) {
                  if (interactive_mode == FALSE) {
                    message("Model averaging enabled.")
                    message(paste0("ans.all$CES: ", ans.all$CES))
                    do.MA <- model_averaging
                  } else {
                    do.MA <- menu(c("no", "yes"), title = "Do you want to calculate the BMD confidence interval by model averaging?") - 
                    1
                  }
                }
                if (do.MA) {
                  if (interactive_mode == FALSE) {
                    ans.all$nr.boot.ma <- num_bootstraps
                  } else {
                    ans.all$nr.boot.ma <- eval(parse(prompt = "\n\ngive number of bootstrap runs for calculating BMD confidence interval based on MA (e.g. 200) >   "))
                  }
                  ans.all$do.MA <- do.MA
                  ans.all$nr.models <- 4
                }
            }
            if (!do.MA && quick.ans < 6) {
              model.options <- c(
                "Exponential model only",
                "Exponential and Hill model",
                "previous option with inverse exponential model added",
                "previous option with lognormal DR model added"
              )
              if (!is.null(model_selection) && model_selection %in% model.options) {
                # If a valid model selection is made, use the corresponding index
                ans.all$nr.models <- match(model_selection, model.options)
              } else {
                model.options <- c("Exponential model only", 
                  "Exponential and Hill model", "previous option with inverse exponential model added", 
                  "previous option with lognormal DR model added")
                cat("Which models do you want to be fitted?\n")
                cat(" ")
                cat(paste(1:4, ":", model.options, "\n"))
                ans.all$nr.models <- eval(parse(prompt = paste(" -------- > ")))
              }
            }
            if (dtype %in% c(5, 15) && !ans.all$no.CI && !ans.all$WAPP && 
                !do.MA && ans.all$CES != 0) 
                ans.all$nruns <- eval(parse(prompt = "\n\ngive the number of bootstrap runs you want to be done \n        for calculating CED CIs for each model (e.g.200) > "))
        }
        if (ans.all$NES.ans == 2) 
            ans.all$CES <- 0.05
        ans.all$no.CI <- F
        if (ans.all$CES == 0) {
            ans.all$CES <- 0.05
            ans.all$no.CI <- T
        }
        data.0.global <- f.remove.NAs(xans = ans.all$xans, yans = Vyans[1], 
            sans = Vsans[1], nans = Vnans[1], covar.no = 0, dfr = ans.all$data.0, 
            output = output)
        if (covar.no > 0) {
            covariate <- data.0.global[, covar.no]
            covariate <- f.remove.blanks(covariate)
            covar.txt <- levels(as.factor(covariate))
        }
        else covar.txt <- character(1)
        CED.matr <- matrix(NA, ncol = (2 * length(covar.txt)))
        dimn.col = character()
        for (ii in 1:length(covar.txt)) dimn.col <- c(dimn.col, 
            paste("BMDL", covar.txt[ii], sep = "-"), paste("BMDU", 
                covar.txt[ii], sep = "-"))
        dimnames(CED.matr)[[2]] <- dimn.col
        ans.all$CED.matr.ma <- matrix(NA, nrow = 0, ncol = 4)
        dimnames(ans.all$CED.matr.ma)[[2]] <- c("endpoint", "subgroup", 
            "BMDL.MA", "BMDU.MA")
        if (ans.all$nr.gr == 1) 
            ans.all$CED.matr.ma <- ans.all$CED.matr.ma[, -2]
        if (ans.all$NES.ans == 2) {
            Mces <- matrix(ncol = 2)
            Msd <- matrix(ncol = 2)
        }
        else {
            Mces <- NA
            Msd <- NA
        }
        ans.all$CI.plt <- FALSE
        ans.all$cont <- f.cont(ans.all)
        if (ans.all$cont) {
            ans.all$plot.type <- 4
            ans.all$plt.mns <- 3
        }
        else ans.all$plot.type <- 2
        CED.all <- list(CED.matr = CED.matr, endpoints = character(), 
            covar.txt = covar.txt, CES = ans.all$CES, Mces = Mces, 
            Msd = Msd, NES.ans = ans.all$NES.ans, nr.models = ans.all$nr.models, 
            quick.ans = quick.ans, PRversion = PRversion, date = date(), 
            Vtrend = logical())
        if (ans.all$NES.ans == 2) 
            CED.all$CES <- NA
        ans.all$model.type <- 1
        ans.all$model.list <- c("EXP", "HILL", "INVEXP", "LOGN")
        CI.matr.ma <- numeric()
        endpoint.ma <- character(0)
        data.0.global <- ans.all$data.0
        first.loop <- T
        ii <- 1
        for (yans.tmp in Vyans) {
            ans.all$yans <- yans.tmp
            if (dtype %in% c(10, 15, 250, 260)) {
                ans.all$sans <- Vsans[ii]
                ans.all$nans <- Vnans[ii]
            }
            ans.all$data.0 <- f.remove.NAs(xans = ans.all$xans, 
                yans = ans.all$yans, sans = ans.all$sans, nans = ans.all$nans, 
                covar.no = 0, dfr = data.0.global, output = output)
            response <- ans.all$data.0[, yans.tmp]
            skip <- FALSE
            if (transf.ans == 3) 
                skip <- f.check.nonneg.num(response, gui, dtype, 
                  quick.ans = 2)
            if (skip) 
                cat("\n\nATTENTION: analysis not possible for response", 
                  ans.all$varnames[yans.tmp], "\n")
            date.tmp <- date()
            if (!skip) {
                if (length(Vyans) > 1) {
                  ans.all$detlim <- 0
                  ans.all$Vdetlim <- 0
                  ans.all$detlim.col <- 0
                  ans.all$covar.no <- covar.no
                }
                ans.all <- f.execute(ans.all, no.plot = T)
                ans.all <- f.clear(ans.all)
                ans.all$twice <- T
                if (first.loop) 
                  ans.all$full.ans <- f.full.ans(ans.all, gui = gui)
                first.loop <- F
                constr.dd <- f.constr.dd(model.ans = 5)
                ans.all$lower.dd <- constr.dd[1]
                ans.all$upper.dd <- constr.dd[2]
                if (!gui && yans.tmp == Vyans[1]) {
                  if (interactive_mode) {
                    if (is.null(lower_dd) && is.null(upper_dd)) {
                        cat(paste("\n\nATTENTION: the constraints on parameter d in the exponential model are set at\n",
                            ans.all$lower.dd, "and", ans.all$upper.dd, "\n"))
                        cat("and adjusted for the other models\n")
                        cat("\ntype 0 if you want to change these constraints,")
                        constr.ans <- as.numeric(readline(prompt = "\notherwise enter any other number > "))
                        if (!is.na(constr.ans) && constr.ans == 0) {
                            ans.all$lower.dd <- as.numeric(readline(prompt = "enter lower bound for parameter d > "))
                            ans.all$upper.dd <- as.numeric(readline(prompt = "enter upper bound for parameter d > "))
                        }
                    } else {
                        if (!is.null(lower_dd)) {
                            ans.all$lower.dd <- lower_dd
                        }
                        if (!is.null(upper_dd)) {
                            ans.all$upper.dd <- upper_dd
                        }
                    }
                  } else {
                    # Non-interactive mode: Use existing defaults for lower_dd and upper_dd if not provided
                    ans.all$lower.dd <- if (is.null(lower_dd)) ans.all$lower.dd else lower_dd
                    ans.all$upper.dd <- if (is.null(upper_dd)) ans.all$upper.dd else upper_dd
                  }
                }
                ans.all$hill <- 0
                ans.all$increase <- 0
                ans.all$Vyans <- Vyans
                ans.all$do.MA <- do.MA
                ans.all <- f.select.con(ans.all, interactive_mode = interactive_mode, results_env = results_env, display_plots = display_plots)
                if (ans.all$quick.ans == 6) {
                  if (!(WAPP || gui)) 
                    f.store.results(ans.all, store.name = 0)
                  return(ans.all)
                }
                do.MA <- ans.all$do.MA
                bmdl.lowest <- numeric()
                bmdu.highest <- numeric()
                if (ans.all$nr.models == 1) {
                  for (qq in (1:nr.gr)) {
                    bmdl.lowest[qq] <- ans.all$EXP$conf.int[qq, 
                      1]
                    bmdu.highest[qq] <- ans.all$EXP$conf.int[qq, 
                      2]
                  }
                }
                if (!ans.all$no.CI && ans.all$nr.models > 0 && 
                  length(ans.all$xans) > 0) {
                  if (ans.all$nr.models > 1) 
                    if (!is.na(ans.all$EXP$conf.int[1, 1]) && 
                      !is.na(ans.all$HILL$conf.int[1, 1])) {
                      nr.gr <- length(ans.all$EXP$conf.int[, 
                        1])
                      nr.gr <- max(nr.gr, length(ans.all$HILL$conf.int[, 
                        1]))
                      if (ans.all$nr.models > 2) 
                        nr.gr <- max(nr.gr, length(ans.all$INVEXP$conf.int[, 
                          1]))
                      if (ans.all$nr.models > 3) 
                        nr.gr <- max(nr.gr, length(ans.all$LOGN$conf.int[, 
                          1]))
                      if (length(ans.all$EXP$conf.int[, 1]) == 
                        1) 
                        ans.all$EXP$conf.int <- matrix(ans.all$EXP$conf.int[1, 
                          ], byrow = TRUE, ncol = 2, nrow = nr.gr)
                      if (length(ans.all$HILL$conf.int[, 1]) == 
                        1) 
                        ans.all$HILL$conf.int <- matrix(ans.all$HILL$conf.int[1, 
                          ], byrow = TRUE, ncol = 2, nrow = nr.gr)
                      for (qq in (1:nr.gr)) {
                        bmdl.lowest[qq] <- min(ans.all$EXP$conf.int[qq, 
                          1], ans.all$HILL$conf.int[qq, 1])
                        bmdu.highest[qq] <- max(ans.all$EXP$conf.int[qq, 
                          2], ans.all$HILL$conf.int[qq, 2])
                      }
                    }
                  if (ans.all$nr.models > 2) 
                    if (!is.na(ans.all$INVEXP$conf.int[1, 1])) {
                      if (length(ans.all$INVEXP$conf.int[, 1]) == 
                        1) 
                        ans.all$INVEXP$conf.int <- matrix(ans.all$INVEXP$conf.int[1, 
                          ], byrow = TRUE, ncol = 2, nrow = nr.gr)
                      for (qq in (1:nr.gr)) {
                        bmdl.lowest[qq] <- min(bmdl.lowest[qq], 
                          ans.all$INVEXP$conf.int[qq, 1])
                        bmdu.highest[qq] <- max(bmdu.highest[qq], 
                          ans.all$INVEXP$conf.int[qq, 2])
                      }
                    }
                  if (ans.all$nr.models > 3) 
                    if (!is.na(ans.all$LOGN$conf.int[1, 1])) {
                      if (length(ans.all$LOGN$conf.int[, 1]) == 
                        1) 
                        ans.all$LOGN$conf.int <- matrix(ans.all$LOGN$conf.int[1, 
                          ], byrow = TRUE, ncol = 2, nrow = nr.gr)
                      for (qq in (1:nr.gr)) {
                        bmdl.lowest[qq] <- min(bmdl.lowest[qq], 
                          ans.all$LOGN$conf.int[qq, 1])
                        bmdu.highest[qq] <- max(bmdu.highest[qq], 
                          ans.all$LOGN$conf.int[qq, 2])
                      }
                    }
                  ans.all$bmdl.lowest <- bmdl.lowest
                  ans.all$bmdu.highest <- bmdu.highest
                  if (ans.all$nr.models > 0) 
                    bmdCI.overall <- rbind(bmdl.lowest, bmdu.highest)
                  bmdCI.overall <- t(bmdCI.overall)
                  if (length(bmdl.lowest) == 1) 
                    covar.txt <- "all"
                  else covar.txt <- dimnames(ans.all$EXP$conf.int)[[1]]
                  bmdCI.overall <- cbind(covar.txt, bmdCI.overall)
                  bmdCI.overall <- as.data.frame(bmdCI.overall)
                  ans.all$bmdCI.overall <- bmdCI.overall
                  if (output) {
                    if (ans.all$NES.ans == 1) 
                      cat("\n\n -----------  CES = ", ans.all$CES, 
                        " --------------------------------- ")
                    if (ans.all$NES.ans == 2 && length(ans.all$bmdCI.overall[, 
                      1]) > 0) 
                      if (output) 
                        cat("\n-----------  endpoint-specific CES", 
                          round(ans.all$CES, 3), "-----------------\n ")
                    if (!do.MA && length(ans.all$bmdCI.overall[, 
                      1]) > 0) {
                      cat("\nThe lowest BMDL and highest BMDU from the fitted models are:\n")
                      print(ans.all$bmdCI.overall)
                      cat(" ---------------------------------------------------\n\n ")
                    }
                  }
                }
                if (length(Vyans) > 1) 
                  if (ans.all$no.CI != T) {
                    CED.all$covar.txt <- ans.all$covar.txt
                    CED.all <- f.CED.all(CED.all, ans.all$y.leg, 
                      ans.all$EXP, ans.all$HILL, ans.all$INVEXP, 
                      ans.all$LOGN, ans.all$TREND)
                  }
                if (interactive_mode == TRUE) {
                  assign("last.fit", ans.all, immed = T, pos = 1)
                }
                if (ans.all$fitted) {
                  if (ans.all$dtype %in% c(1, 5)) 
                    ans.all$plt.mns <- 1
                  ans.all <- f.plot.gui(ans.all, display_plots = display_plots, results_env = results_env)
                  if (!WAPP && length(ans.all$xans) > 1 && ans.all$nr.aa > 
                    1) {
                    f.press.key.to.continue()
                    f.plot.all(ans.all, sep = T, display_plots = display_plots)
                  }
                }
                if (length(ans.all$xans) > 1) {
                  ans.all$conf.int <- f.confint.doseaddition(ans.all)
                  cat("\n the confidence intervals (per model) are:\n")
                  print(ans.all$conf.int)
                }
                if (do.MA) {
                  if (ans.all$TREND) {
                    ans.all$MA.running <- TRUE
                    cat("\n\nCalculating confidence intervals by model averaging, this may make some time ....\n\n")
                    ans.all <- f.boot.ma(ans.all, interactive_mode = interactive_mode, results_env = results_env, display_plots = display_plots)
                    cat("\nThe model-average BMD confidence interval is:\n")
                    print(ans.all$MA$conf.int.ma)
                    cat("\n")
                    CI.ma <- ans.all$MA$conf.int.ma
                    CI.ma <- cbind(rep(ans.all$y.leg, ans.all$nr.gr), 
                      CI.ma)
                    dimnames(CI.ma)[[2]][1] <- "endpoint"
                    ans.all$CED.matr.ma <- rbind(ans.all$CED.matr.ma, 
                      CI.ma)
                    CI.matr.ma <- rbind(CI.matr.ma, ans.all$MA$CI.row.ma)
                    if (interactive_mode == FALSE) {
                      if (length(ans.all$Vyans) > 1) {
                        message("Vyans is greater than 1")
                        ans.all$res.name <- ans.all$varname[ans.all$yans]
                        assign(paste(ans.all$res.name, "model_averaging"), ans.all, envir = results_env)
                      } else {
                        ans.all$res.name <- ans.all$varname[ans.all$yans]
                        results_env$model_averaging <- ans.all
                      } 
                    }
                  } else {
                    ans.all$MA <- list()
                    CI.matr.ma <- rbind(CI.matr.ma, rep(NA, 2 * 
                      ans.all$nr.gr))
                    ans.all$MA$conf.int.ma <- CI.matr.ma
                    cat("\n\nNo Model Averaging applied, due to nonsignificant trend in the data\n\n")
                  }
                  endpoint.ma <- c(endpoint.ma, ans.all$y.leg)
                }
                if (length(ans.all$HILL) == 0) 
                  ans.all <- f.move.sublist(ans.all, ans.all$EXP)
                if (!gui && !WAPP) {
                  if (length(ans.all$notes) > 1) {
                    print("f.quick.con")
                    cat(ans.all$notes)
                    f.press.key.to.continue()
                  }
                  if (!(WAPP || gui) && length(Vyans) > 1 && 
                    interrupt) 
                    store.name <- f.store.results(ans.all, store.name = 0)
                }
            }
            if (!gui) 
                cat(paste("\n\nend of analysis for response: ", 
                  ans.all$y.leg))
            if (gui) 
                ans.all$date <- date()
            ans.all$MA.running <- FALSE
            ii <- ii + 1
        }
        if (ans.all$do.MA) {
            CED.matr.ma.plt <- data.frame(endpoint = endpoint.ma, 
                CI = CI.matr.ma)
            ans.all$CED.matr.ma.plt <- CED.matr.ma.plt
        }
        if (ans.all$no.CI == F && length(Vyans) == 1 && length(ans.all$EXP$conf.int[, 
            1]) > 1 && length(ans.all$xans) == 1) 
            f.plot.CI(ans.all, display_plots = display_plots)
        if (ans.all$no.CI == F & length(Vyans) > 1) {
            if (sum(!is.na(CED.all$CED.matr[, 1]) > 0)) {
                CED.all$CED.matr <- CED.all$CED.matr[-1, ]
                dimnames(CED.all$CED.matr)[[1]] <- CED.all$endpoints
                dimnames(CED.all$Mces)[[1]] <- NULL
                CED.all$Mces <- as.data.frame(CED.all$Mces[-1, 
                  ])
                dimnames(CED.all$Msd)[[1]] <- NULL
                CED.all$Msd <- as.data.frame(CED.all$Msd[-1, 
                  ])
                CED.all$max.x <- max(ans.all$x)
                lst.name <- 0
                if (WAPP) 
                  graphics.off()
                if (output) {
                  cat("\ncalculations started at:", date.tmp)
                  cat("\nand ended at:", date())
                  cat("\n")
                }
                ans.all$data.0 <- data.0.global
                ans.all$means.plt.lst <- NULL
                ans.all$Vyans <- Vyans
                ans.all$CED.all <- CED.all
                f.plot.CED(ans.all, WAPP = WAPP, plotprefix = ans.all$plotprefix, 
                  svg.plots = ans.all$svg.plots, display_plots = display_plots, interactive_mode = interactive_mode)
            }
            else cat("\n\nNone of the endpoints showed a significant trend\n")
        }
        if (!gui) 
            ans.all$NES.ans <- 1
        if (exists("track")) 
            print("f.quick.con:  END")
        return(ans.all)
  })
}



f.lik.con <- function(theta, x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
    sd2.log, nn, Vdetlim, CES, twice = T, ttt = 0, trace.tmp = F, 
    fct4 = 1, fct5 = 1, cens.up = NA, lb = -Inf, ub = Inf, par.tmp, 
    increase = increase, x.mn = NA, ref.lev = ref.lev, ans.m6.sd = 1, 
    sign.q, Mx = 0, x1 = 0, x2 = 0, cc.inf, incr.decr.no = 0) {
    if (model.ans == -42) {
        print("f.lik.con:  begin")
        print(model.ans)
        print(theta)
    }
    if (sum(is.na(theta)) > 0) {
        if (trace.tmp) {
            cat("\nnote: problem in f.lik.con: NAs in theta\n")
            print(theta)
        }
        return(NA)
    }
    if (0) 
        if ((length(fct3) > 1) & (length(fct3) != length(x))) {
            print("f.lik.con")
            cat("\nATTENTION: fct3 incorrect\n\n")
            print(length(fct3))
            stop()
        }
    variance <- 0
    for (jj in 1:max(fct3)) variance <- variance + theta[jj] * 
        (fct3 == jj)
    if (!model.ans %in% c(6, 47) || ans.m6.sd != 2) 
        regr.par <- theta[(max(fct3) + 1):length(theta)]
    else regr.par <- theta
    if (model.ans == 11) {
        yy <- y
        yy[y == -1000] <- Vdetlim[y == -1000]
        yy[y == -2000] <- cens.up
    }
    if (sum(!is.finite(theta)) > 0) {
        theta[!is.finite(theta)] <- par.tmp[!is.finite(theta)]
    }
    if (length(Mx) != 1) {
        nr.dosecol <- length(Mx[1, ])
        RPF.vec <- regr.par[(length(regr.par) - nr.dosecol + 
            2):length(regr.par)]
        x <- Mx %*% c(1, RPF.vec)
    }
    if (sum(theta <= lb) > 0) {
        theta[theta < lb] <- 1.1 * par.tmp[theta < lb]
    }
    if (sum(theta >= ub) > 0) {
        theta[theta > ub] <- 0.9 * par.tmp[theta > ub]
    }
    if (dtype %in% c(1, 5, 25, 26)) 
        expect <- f.expect.con(model.ans, x, regr.par, fct1 = fct1, 
            fct2 = fct2, fct3 = fct3, fct4 = fct4, fct5 = fct5, 
            CES = CES, twice = twice, ttt = ttt, y = yy, trace.expect = F, 
            increase = increase, x.mn = x.mn, ref.lev = ref.lev, 
            sign.q = sign.q, ans.m6.sd = ans.m6.sd, x1 = x1, 
            x2 = x2, cc.inf = cc.inf, incr.decr.no = incr.decr.no)
    if (dtype == 1 | dtype == 5) {
        if (trace.tmp & sum(is.na(expect)) > 0) {
            cat("\nnote: NAs in predicted response\nat parameter values:\n")
            cat(signif(theta, 4), "\n")
        }
        expect <- logb(expect)
        y.log <- logb(y)
        y.log[y < 0] <- 0
        score1 <- (y > 0) * (-0.5 * logb(2 * pi * variance) - 
            ((y.log - expect)^2)/(2 * variance))
        score.detlim <- logb(pnorm((logb(Vdetlim) - expect)/sqrt(variance)))
        score.detlim[!is.finite(score.detlim)] <- 0
        score.censup <- logb(1 - pnorm((logb(cens.up) - expect)/sqrt(variance)))
        score.censup[!is.finite(score.censup)] <- 0
        score2 <- (y == -1000) * score.detlim + (y == -2000) * 
            score.censup
        score <- score1 + score2
    }
    if (dtype == 25) {
        score1 <- (y > 0) * (-0.5 * logb(2 * pi * variance) - 
            ((y - expect)^2)/(2 * variance))
        score.detlim <- logb(pnorm(((Vdetlim) - expect)/sqrt(variance)))
        score.detlim[!is.finite(score.detlim)] <- 0
        score.censup <- logb(1 - pnorm(((cens.up) - expect)/sqrt(variance)))
        score.censup[!is.finite(score.censup)] <- 0
        score2 <- (y == -1000) * score.detlim + (y == -2000) * 
            score.censup
        score <- score1 + score2
    }
    if (dtype == 26) {
        expect <- sqrt(expect)
        y.sqrt <- sqrt(y)
        score1 <- (y > 0) * (-0.5 * logb(2 * pi * variance) - 
            ((y.sqrt - expect)^2)/(2 * variance))
        score.detlim <- logb(pnorm((sqrt(Vdetlim) - expect)/sqrt(variance)))
        score.detlim[!is.finite(score.detlim)] <- 0
        score.censup <- logb(1 - pnorm((sqrt(cens.up) - expect)/sqrt(variance)))
        score.censup[!is.finite(score.censup)] <- 0
        score2 <- (y == -1000) * score.detlim + (y == -2000) * 
            score.censup
        score <- score1 + score2
    }
    if (dtype %in% c(10, 15, 250, 260)) 
        expect <- f.expect.con(model.ans, x, regr.par, fct1 = fct1, 
            fct2 = fct2, fct3 = fct3, fct4 = fct4, fct5 = fct5, 
            CES = CES, twice = twice, ttt = ttt, y = y, increase = increase, 
            x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
            ans.m6.sd = ans.m6.sd, x1 = x1, x2 = x2, cc.inf = cc.inf, 
            incr.decr.no = incr.decr.no)
    if (dtype %in% c(10, 15)) {
        if (model.ans != 11) {
            expect <- logb(expect)
            dum <- nn * (mn.log - expect)^2 + (nn - 1) * sd2.log
        }
        else dum <- (nn - 1) * sd2.log
        score <- -(nn * logb(sqrt(2 * pi * variance)) + dum/(2 * 
            variance))
    }
    if (dtype == 250) {
        if (model.ans != 11) {
            dum <- nn * (mn.log - expect)^2 + (nn - 1) * sd2.log
        }
        else dum <- (nn - 1) * sd2.log
        score <- -(nn * logb(sqrt(2 * pi * variance)) + dum/(2 * 
            variance))
    }
    if (dtype == 260) {
        if (model.ans != 11) {
            expect <- sqrt(expect)
            dum <- nn * (mn.log - expect)^2 + (nn - 1) * sd2.log
        }
        else dum <- (nn - 1) * sd2.log
        score <- -(nn * logb(sqrt(2 * pi * variance)) + dum/(2 * 
            variance))
    }
    if (model.ans == -51) {
        print("f.lik.con END")
        print(theta)
        print(sum(score))
    }
    return(-sum(score))
}



f.expect.con <- function(model.ans, x, regr.par = 0, fct1 = 1, fct2 = 1, fct3 = 1, 
    fct4 = 1, fct5 = 1, name = F, CES = NA, twice = T, ttt = 0, 
    yy = 0, trace.expect = F, increase, x.mn = NA, ref.lev, ans.m6.sd = 1, 
    par.start = NA, sign.q = 0, x1, x2, opposing = 1, cc.inf, 
    incr.decr.no = 0, cont = TRUE) {
    if (exists("track2")) 
        print("f.expect.con")
    if (0) 
        if (name == !T) {
            print("f.expect.con ---    ")
            print(model.ans)
            print(x)
        }
    if (name) {
        model.names <- c("E1: y = a", "E2: y = a*exp(bx)", "E3: y = a*exp(bx^d)", 
            "E4: y = a*[c-(c-1)exp(-bx)]", "E5: y = a*[c-(c-1)exp(-bx^d)]", 
            "E5b: y = a*(c^(1 - exp(-(x/b)^d)))", "E2 y = a*exp(x/b)", 
            "E3 y = a*exp((x/b)^d)", "E4 y = a * [c - (c-1)exp(-(x/b))]", 
            "E5 y = a * [c - (c-1)exp(-(x/b)^d)]", "Full model: y = group mean", 
            "E2-CED: y = a*exp(bx)", "E3-CED: y = a*(c^(1 - exp(-(x/b)^d))) with c='Inf'", 
            "E4-CED: y = a * [c-(c-1)exp(-bx)]", "E5-CED: y = a*(c^(1 - exp(-(x/b)^d)))", 
            "BMDratio model (M10)", "H2: a * (1 - x/(b+x))", 
            "H3: a * (1 - x^d/(b^d+x^d))", "H4: a * (1 + (c-1)x/(b+x))", 
            "H5: a * (1 + (c-1)x^d/(b^d+x^d))", "H5b: y = a*(c^(x^d/(b^d + x^d)))", 
            "H2-CED: a * (1 - x/(b+x))", "H3-CED: a * c ^(x^d/(b^d+x^d)) with c='Inf'", 
            "H4-CED: a * (1 + (c-1)x/(b+x))", "H5-CED: a * c ^(x^d/(b^d+x^d))", 
            "y = c + a*x^b", "y = c + exp((x/b)^a)", "y = c + b(x-a)^d*(x>a)", 
            "y = a + c * (1 - exp( - (x/b)^d))", "y = a + b*x^c", 
            "y = a {c1^(1 - exp(-(x/b1)^d)) * c2^(1 - exp(-(x/b2)^d))}", 
            "y = a {[cc-(cc-1)exp(-(x/b1)^d1)] * [cc - (cc-1)exp(-(x/b2)^d2)]}", 
            "y = a {c1^(1 - exp(-(x/b1)^d1)) * c2^(1 - exp(-(x/b2)^d2))}", 
            "y = a + bx + cx^2 + dx^3", "y = a / b^x", "y = a * b^x", 
            "y =  b*exp(c(x-a)^d) * (x>a) + b * (x<a)", "select model 5 from nested families of models", 
            "select model 3 or model 5 from nested families of models", 
            "select model 3 from nested families of models", 
            "y = ax^b", "effect multiplication based on model 6", 
            "hill-additive: y = a + cx^d/(b^d+x^d)", "E2 with independent a and b: y = exp(a + b(x - mean.x))", 
            "y = a {[c1-(c1-1)exp(-(x/b)^d)] * [c2 - (c2-1)exp(-(x/b)^d)]}", 
            "E5-CED in terms of RPFs", "E5b in terms of CED, with c replaced by qs according to ES theory", 
            "CxT model 6 with b = aa.t * exp( - (time/bb.t)^dd.t)", 
            "E5b-ED50: y = a*(c^(1 - exp(-ln(2)(x/ED50)^d)))", 
            "CxT model 15 with CED = aa.t * exp( - (time/bb.t)^dd.t)", 
            "inverse exponential in terms of CED with c='Inf'", 
            "inverse exponential in terms of CED", "LN DR model in terms of CED with c='Inf'", 
            "LN DR model in terms of CED", "AUC-external dose: y = x * (Km + 0.5*x/V) / Vmax", 
            "SS-external dose: y = Km * x/V / (Vmax - x/V)", 
            "y = a {cc^(1 - exp(-(x/b1)^d1)) * exp(-(x/b2)^d2)}", 
            "y = (FD / V) * ka / (ka - ke) * [exp(- ke*t) - exp(- ka*t)]", 
            "y = C1 * exp(- lab1*t) + C2 * exp(-lab2*t) - (C1 + C2) * exp(- ka*t)", 
            "logE5-canon: y = exp[ a(1 + (c - 1)(1 - exp(-bx^d)) ]", 
            "logE5-noncan: y = exp[ a(1 + (c - 1)(1 - exp(-bx^d)) ]", 
            "EFSA probit model, noncanon", "EFSA probit model y = exp(a * pnorm(c + b * x^d))", 
            "two comp. PBPK lipoph", "simple two-comp model")
        return(model.names)
    }
    y.expect <- nr.aa <- max(fct1)
    nr.bb <- max(fct2)
    nr.var <- max(fct3)
    nr.cc <- max(fct4)
    nr.dd <- max(fct5)
    if (model.ans %in% c(31, 33, 57)) {
        nr.b1 <- 1
        nr.b2 <- nr.bb
    }
    if (model.ans == 57) 
        nr.b1 <- nr.dd
    if (model.ans == 47 || (model.ans == 6 && ans.m6.sd == 2)) {
        sd.tmp <- sqrt(regr.par[1:nr.var])
        sd0 <- rep(0, length(x))
        for (ii in (1:nr.var)) sd0 <- sd0 + sd.tmp[ii] * (fct3 == 
            ii)
        regr.par <- regr.par[-(1:nr.var)]
        qq <- regr.par[nr.aa + nr.bb + 1]
        sign.q0 <- rep(0, length(x))
        for (ii in (1:nr.var)) sign.q0 <- sign.q0 + sign.q[ii] * 
            (fct3 == ii)
    }
    if (model.ans != 11) {
        nrp <- length(regr.par)
        aa0 <- rep(0, length(x))
        aa.tmp <- regr.par[1:nr.aa]
        if (nr.aa == 1) 
            aa0 <- aa0 + aa.tmp[1]
        else for (ii in (1:nr.aa)) aa0 <- aa0 + aa.tmp[ii] * 
            (fct1 == ii)
        bb0 <- rep(0, length(x))
        bb.tmp <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
        if (model.ans %in% c(31, 33, 57)) 
            bb.tmp <- regr.par[(nr.aa + nr.b1 + 1):(nr.aa + nr.b1 + 
                nr.b2)]
        for (jj in (1:nr.bb)) bb0 <- bb0 + bb.tmp[jj] * (fct2 == 
            jj)
        if (model.ans == 57) {
            bb2.0 <- bb0
            bb0 <- rep(0, length(x))
            bb.tmp <- regr.par[(nr.aa + 1):(nr.aa + nr.b1)]
            for (jj in (1:nr.b1)) bb0 <- bb0 + bb.tmp[jj] * (fct5 == 
                jj)
            bb1.0 <- bb0
        }
        par3 <- regr.par[nr.aa + nr.bb + 1]
        if (length(par3) == 0 || is.na(par3)) 
            par3 <- 0
        par4 <- regr.par[nr.aa + nr.bb + nr.cc + 1]
        cc0 <- par3
        dd0 <- par4
        if (model.ans %in% c(3, 8, 13, 18, 23, 51, 53, 65)) {
            dd0 <- par3
            if (max(fct5) > 1) {
                dd0 <- rep(0, length(x))
                dd.tmp <- regr.par[(nr.aa + nr.bb + 1):length(regr.par)]
                for (kk in (1:nr.dd)) dd0 <- dd0 + dd.tmp[kk] * 
                  (fct5 == kk)
            }
        }
        if (model.ans %in% c(5, 6, 10, 15, 20, 21, 25, 41, 42, 
            46, 52, 54)) {
            dd0 <- par4
            if (max(fct5) > 1) {
                dd0 <- rep(0, length(x))
                dd.tmp <- regr.par[(nr.aa + nr.bb + nr.cc + 1):length(regr.par)]
                for (kk in (1:nr.dd)) dd0 <- dd0 + dd.tmp[kk] * 
                  (fct5 == kk)
            }
        }
        if (model.ans %in% c(4, 5, 6, 9, 10, 14, 15, 19, 20, 
            24, 25, 41, 42, 46, 52, 54, 55, 56)) {
            if (max(fct4) > 1) {
                cc0 <- rep(0, length(x))
                cc.tmp <- regr.par[(nr.aa + nr.bb + 1):(nr.aa + 
                  nr.bb + nr.cc)]
                for (kk in (1:nr.cc)) cc0 <- cc0 + cc.tmp[kk] * 
                  (fct4 == kk)
            }
        }
        if (model.ans %in% c(31, 33, 57)) {
            cc0 <- rep(0, length(x))
            cc.tmp <- regr.par[(nr.aa + nr.b1 + nr.b2 + 1):(nr.aa + 
                nr.b1 + nr.b2 + nr.cc)]
            for (kk in (1:nr.cc)) cc0 <- cc0 + cc.tmp[kk] * (fct4 == 
                kk)
        }
    }
    if (model.ans %in% c(49, 50)) {
        gg <- regr.par[nr.aa + nr.bb + nr.cc + nr.dd + 1]
        hh <- regr.par[nr.aa + nr.bb + nr.cc + nr.dd + 2]
    }
    switch(model.ans, y.expect <- aa0, y.expect <- aa0 * exp(bb0 * 
        x + cc0 * ttt), y.expect <- aa0 * exp(bb0 * (x^dd0)), 
        y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-bb0 * x)), 
        y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-bb0 * (x^dd0))), 
        {
            if (ans.m6.sd == 1) y.expect <- exp(log(aa0) + log(cc0) * 
                (1 - exp(-(x/bb0)^dd0)))
            if (ans.m6.sd == 2) y.expect <- exp(log(aa0) + sign.q0 * 
                qq * sd0 * (1 - exp(-(x/bb0)^dd0)))
        }, y.expect <- aa0 * exp(x/bb0), {
            if (bb0[1] > 0) y.expect <- aa0 * exp((x/bb0)^dd0) else y.expect <- aa0 * 
                exp(-(x/-bb0)^dd0)
        }, y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-x/bb0)), 
        y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-(x/bb0)^dd0)), 
        {
            x.gr <- levels(factor(x))
            x.fact <- factor(x)
            y.tmp <- rep(NA, length(x))
            y.expect <- rep(0, length(x))
            if (nr.var > 1 & nr.aa == 1 & nr.bb == 1) {
                for (jj in (1:nr.var)) for (ii in (1:length(x.gr))) {
                  y.tmp <- yy[x.fact == x.gr[ii] & fct3 == jj]
                  if (length(y.tmp) > 0) {
                    y.mn <- exp(mean(logb(y.tmp)))
                    y.expect <- y.expect + y.mn * (x.fact == 
                      x.gr[ii]) * (fct3 == jj)
                  }
                  y.mn <- y.tmp
                }
            } else if (twice) for (jj in (1:nr.bb)) for (ii in (1:length(x.gr))) {
                y.tmp <- yy[x.fact == x.gr[ii] & fct2 == jj]
                if (length(y.tmp) > 0) {
                  y.mn <- exp(mean(logb(y.tmp)))
                  y.expect <- y.expect + y.mn * (x.fact == x.gr[ii]) * 
                    (fct2 == jj)
                }
            } else if (!twice) for (jj in (1:nr.aa)) for (kk in (1:nr.bb)) for (ii in (1:length(x.gr))) {
                y.tmp <- yy[x.fact == x.gr[ii] & fct1 == jj & 
                  fct2 == kk]
                if (length(y.tmp) > 0) {
                  y.mn <- exp(mean(logb(y.tmp)))
                  y.expect <- y.expect + y.mn * (x.fact == x.gr[ii]) * 
                    (fct1 == jj) * (fct2 == kk)
                }
            }
        }, {
            y.expect <- aa0 * (CES + 1)^(x/bb0)
        }, {
            if (increase == 1) cc0 <- cc.inf
            if (increase == -1) cc0 <- 1/cc.inf
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, cont = cont)
            y.expect <- log(aa0) + log(cc0) * (1 - exp(-(x/dum)^dd0))
            y.expect <- exp(y.expect)
        }, y.expect <- aa0 * (cc0 - (cc0 - 1) * ((CES + 1 - cc0)/(1 - 
            cc0))^(x/bb0)), {
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, incr.decr.no = incr.decr.no, cont = cont)
            y.expect <- log(aa0) + log(cc0) * (1 - exp(-(x/dum)^dd0))
            y.expect <- exp(y.expect)
        }, {
            if (!exists("CES.16")) {
                cat("\n\nATTENTION:  CES.16 not defined; CES2 is set to 0.10 and CES1 to 0.05\n\n")
                CES2 <- 0.1
                CES1 <- 0.05
                f.press.key.to.continue()
            } else {
                CES2 <- CES.16$CES2
                CES1 <- CES.16$CES1
            }
            if (par4 <= 1) print(c("f.expect.con, value of BMD-ratio:", 
                par4))
            dd <- f.uniroot.BMDratio(ratio = par4, cc = cc0, 
                CES1 = CES1, CES2 = CES2)
            y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-(x/bb0)^dd))
            if (is.na(dd)) y.expect <- rep(0, length(x))
        }, {
            y.expect <- aa0 * (1 - x/(bb0 + x))
        }, {
            if (increase == 1) y.expect <- aa0 * (1 - x^dd0/(-(-bb0)^dd0 + 
                x^dd0))
            if (increase == -1) y.expect <- aa0 * (1 - x^dd0/((bb0)^dd0 + 
                x^dd0))
        }, {
            y.expect <- aa0 * (1 + ((cc0 - 1) * x)/(bb0 + x))
        }, {
            y.expect <- aa0 * (1 + (cc0 - 1) * x^dd0/(bb0^dd0 + 
                x^dd0))
        }, {
            y.expect <- log(aa0) + log(cc0) * (x^dd0/(bb0^dd0 + 
                x^dd0))
            y.expect <- exp(y.expect)
        }, {
            dum <- f.bb.con(model.ans, cc = NA, dd = NA, CED = bb0, 
                CES, cont = cont)
            y.expect <- aa0 * (1 - x/(sign(dum) * (abs(dum) + 
                x)))
            y.expect <- aa0 * (1 - x/(dum + x))
        }, {
            if (increase == 1) cc0 <- cc.inf
            if (increase == -1) cc0 <- 1/cc.inf
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, cont = cont)
            y.expect <- log(aa0) + log(cc0) * (x^dd0/(dum^dd0 + 
                x^dd0))
            y.expect <- exp(y.expect)
        }, {
            dum <- f.bb.con(model.ans, cc = cc0, dd = NA, CED = bb0, 
                CES, cont = cont)
            y.expect <- aa0 * (1 + ((cc0 - 1) * x)/(dum + x))
        }, {
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, incr.decr.no = incr.decr.no, cont = cont)
            y.expect <- log(aa0) + log(cc0) * (x^dd0/(dum^dd0 + 
                x^dd0))
            y.expect <- exp(y.expect)
        }, {
            y.expect <- cc0 + aa0 * x^bb0
        }, y.expect <- cc0 + exp(x/bb0)^aa0, y.expect <- cc0 + 
            (bb0 * abs((x - aa0))^dd0) * (x >= aa0), y.expect <- aa0 + 
            cc0 * (1 - exp(-(x/bb0)^dd0)), {
            y.expect <- aa0 + bb0 * x^cc0
        }, {
            b1 <- regr.par[nr.aa + 1]
            c2 <- regr.par[nr.aa + 1 + nr.b2 + nr.cc + 1]
            dd <- regr.par[nr.aa + 1 + nr.b2 + nr.cc + 2]
            if (0) {
                print(" f.expect.con")
                print(regr.par)
                print(aa0[1:5])
                print(b1[1:5])
                print(bb0[1:20])
                print(cc0[1:20])
                print(c2[1:5])
                print(dd[1:5])
            }
            y.expect.1 <- log(cc0) * (1 - exp(-(x/b1)^dd))
            y.expect.2 <- log(c2) * (1 - exp(-(x/bb0)^dd))
            y.expect <- exp(log(aa0) + y.expect.1 + y.expect.2)
        }, {
            b2 <- regr.par[nr.aa + nr.bb + 1]
            cc <- regr.par[nr.aa + nr.bb + 2]
            d1 <- regr.par[nr.aa + nr.bb + 3]
            d2 <- regr.par[nr.aa + nr.bb + 4]
            y.expect.1 <- cc - (cc - 1) * exp(-(x/bb0)^d1)
            y.expect.2 <- cc - (cc - 1) * exp(-(x/b2)^d1)
            y.expect <- aa0 * (y.expect.1 * y.expect.2)
        }, {
            b1 <- regr.par[nr.aa + 1]
            c2 <- regr.par[nr.aa + 1 + nr.b2 + nr.cc + 1]
            d1 <- regr.par[nr.aa + 1 + nr.b2 + nr.cc + 2]
            d2 <- regr.par[nr.aa + 1 + nr.b2 + nr.cc + 3]
            if (0) {
                print(" f.expect.con")
                print(regr.par)
                print(aa0[1:5])
                print(b1[1:5])
                print(bb0[1:20])
                print(cc0[1:20])
                print(c2[1:5])
                print(d1[1:5])
                print(d2[1:5])
            }
            y.expect.1 <- log(cc0) * (1 - exp(-(x/b1)^d1))
            y.expect.2 <- log(c2) * (1 - exp(-(x/bb0)^d2))
            y.expect <- exp(log(aa0) + y.expect.1 + y.expect.2)
        }, {
            y.expect <- aa0 + bb0 * x + cc0 * x^2 + dd0 * x^3
        }, y.expect <- aa0/bb0^x, y.expect <- aa0 * bb0^x, {
            tmp1 <- bb0 * (x < aa0) + 0 * (x >= aa0)
            tmp2 <- bb0 * exp(cc0 * (x - aa0)^dd0) * (x >= aa0)
            tmp2[is.na(tmp2)] <- 0
            y.expect <- tmp1 + tmp2
        }, cat(""), cat(""), cat(""), y.expect <- aa0 * x^bb0, 
        {
            aa <- regr.par[1]
            bb1 <- regr.par[2]
            bb2 <- regr.par[3]
            cc1 <- regr.par[4]
            cc2 <- regr.par[5]
            dd <- regr.par[6]
            y.expect.1 <- log(aa) + log(cc1) * (1 - exp(-(x1/bb1)^dd)) * 
                (x2 == 0)
            y.expect.2 <- log(aa) + log(cc1) * (1 - exp(-(x1/bb1)^dd)) * 
                (x1 == 0)
            y.expect.bis <- y.expect.1 * (x1 != 0) * y.expect.2 * 
                (x2 != 0)
            y.expect <- exp(y.expect.1 + y.expect.2 + y.expect.bis)
        }, {
            y.expect <- aa0 + cc0 * x^dd0/(bb0^dd0 + x^dd0)
        }, {
            y.expect <- exp(aa0 + bb0 * (x - x.mn))
        }, {
            c1 <- regr.par[nr.aa + nr.bb + 1]
            c2 <- regr.par[nr.aa + nr.bb + 2]
            dd <- regr.par[nr.aa + nr.bb + 3]
            y.expect.1 <- c1 - (c1 - 1) * exp(-(x/bb0)^dd)
            y.expect.2 <- c2 - (c2 - 1) * exp(-(x/bb0)^dd)
            y.expect <- aa0 * (y.expect.1 * y.expect.2)
        }, {
            CED.pars <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
            CED.ref <- CED.pars[ref.lev]
            bb0 <- CED.ref/bb0
            bb0[fct2 == ref.lev] <- CED.ref
            dum <- ((x/bb0)^dd0) * (-log(1 - log(CES + 1)/log(cc0)))
            y.expect <- log(aa0) + log(cc0) * (1 - exp(-dum))
            y.expect <- exp(y.expect)
        }, {
            y.expect <- exp(log(aa0) + sign.q0 * qq * sd0 * (1 - 
                exp(x^dd0 * log(1 - CES)/bb0^dd0)))
        }, {
            aa0 <- regr.par[1]
            cc0 <- regr.par[2]
            dd0 <- regr.par[3]
            aa.t <- regr.par[4]
            bb.t <- regr.par[5]
            dd.t <- regr.par[6]
            bb0 <- aa.t * exp(-(ttt/bb.t)^dd.t)
            y.expect <- exp(log(aa0) + log(cc0) * (1 - exp(-(x/bb0)^dd0)))
        }, {
            y.expect <- exp(log(aa0) + log(cc0) * (1 - exp(-log(2) * 
                (x/bb0)^dd0)))
        }, {
            aa0 <- regr.par[1]
            cc0 <- regr.par[2]
            dd0 <- regr.par[3]
            aa.t <- regr.par[4]
            bb.t <- regr.par[5]
            dd.t <- regr.par[6]
            CED0 <- aa.t * exp(-(ttt/bb.t)^dd.t)
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = CED0, 
                CES, cont = cont)
            y.expect <- exp(log(aa0) + log(cc0) * (1 - exp(-(x/dum)^dd0)))
        }, {
            if (increase == 1) cc0 <- cc.inf
            if (increase == -1) cc0 <- 1/cc.inf
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, cont = cont)
            y.expect <- log(aa0) + log(cc0) * exp(-(x/dum)^(-dd0))
            y.expect <- exp(y.expect)
        }, {
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, incr.decr.no = incr.decr.no, cont = cont)
            y.expect <- log(aa0) + log(cc0) * exp(-(x/dum)^(-dd0))
            y.expect <- exp(y.expect)
        }, {
            if (increase == 1) cc0 <- cc.inf
            if (increase == -1) cc0 <- 1/cc.inf
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, cont = cont)
            y.expect <- log(aa0) + log(cc0) * pnorm(log(dum) + 
                dd0 * log(x))
            y.expect <- exp(y.expect)
        }, {
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, incr.decr.no = incr.decr.no, cont = cont)
            y.expect <- log(aa0) + log(cc0) * pnorm(log(dum) + 
                dd0 * log(x))
            y.expect <- exp(y.expect)
        }, {
            y.expect <- x * (aa0 + 0.5 * x/cc0)/bb0
        }, {
            y.expect <- aa0 * x/cc0/(bb0 - x/cc0)
        }, {
            d1 <- regr.par[nr.aa + nr.b1 + nr.b2 + nr.cc + 1]
            d2 <- regr.par[nr.aa + nr.b1 + nr.b2 + nr.cc + 2]
            y.expect.1 <- log(cc0) * (1 - exp(-(x/bb1.0)^d1))
            y.expect.2 <- log(exp(-(x/bb2.0)^d2))
            y.expect <- exp(log(aa0) + y.expect.1 + y.expect.2)
        }, {
            FD.V <- regr.par[1]
            ka <- regr.par[2]
            ke <- regr.par[3]
            y.expect <- FD.V * (ka/(ka - ke)) * (exp(-ke * x) - 
                exp(-ka * x))
        }, {
            C1 <- regr.par[1]
            C2 <- regr.par[2]
            ka <- regr.par[3]
            lab1 <- regr.par[4]
            lab2 <- regr.par[5]
            y.expect <- C1 * exp(-lab1 * x) + C2 * exp(-lab2 * 
                x) - (C1 + C2) * exp(-ka * x)
        }, {
            aa.can <- aa0
            cc.can <- cc0
            aa0 <- log(aa.can)
            cc0 <- 1 + log(cc.can)/log(aa.can)
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, incr.decr.no = incr.decr.no, cont = cont, 
                aa = aa0)
            y.expect <- exp(aa0 * (1 + (cc0 - 1) * (1 - exp(-dum * 
                (x^dd0)))))
        }, {
            dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                CES, incr.decr.no = incr.decr.no, cont = cont, 
                aa = aa0)
            y.expect <- exp(aa0 * (1 + (cc0 - 1) * (1 - exp(-dum * 
                (x^dd0)))))
        }, {
            y.expect <- exp(aa0 * pnorm(cc0 + bb0 * x^dd0))
        }, {
            aa.can <- aa0
            cc.can <- cc0
            aa0 <- log(aa.can * cc.can)
            cc0 <- qnorm(log(aa.can)/log(aa.can * cc.can))
            dum <- f.bb.con(model.ans = 63, cc = cc0, dd = dd0, 
                CED = bb0, CES, incr.decr.no = incr.decr.no, 
                cont = cont, aa = aa0)
            y.expect <- exp(aa0 * pnorm(cc0 + dum * x^dd0))
        }, {
            D <- regr.par[1]
            V <- regr.par[2]
            lab.min <- regr.par[3]
            lab.plus <- regr.par[4]
            qk <- regr.par[5]
            dum = (D/V)/(lab.min - lab.plus)
            y.expect = dum * ((lab.min + qk) * exp(lab.plus * 
                x) - (lab.plus + qk) * exp(lab.min * x))
        }, {
            cc <- aa0
            k1 <- bb0
            ke <- dd0
            dum = cc * k1/(k1 - ke)
            y.expect = dum * (exp(-ke * x) - exp(-k1 * x))
        })
    if (exists("track2")) 
        print("f.expect.con:  END")
    return(y.expect)
}



f.constr.con <- function(ans.all, tmp.quick = F) {
    if (exists("track")) 
        print("f.constr.con")
    ans.all$tmp.quick <- tmp.quick
    with(ans.all, {
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        nr.cc <- max(fct4)
        nr.dd <- max(fct5)
        if (dtype != 6) 
            nr.var <- max(fct3)
        lower.var <- 1e-06
        upper.var <- 10
        if (dtype %in% c(25, 250)) 
            upper.var <- 1e+10
        if (dtype %in% c(26, 260)) 
            upper.var <- 10000
        lower.RPF <- 1e-12
        upper.RPF <- 1e+12
        if (tmp.quick && model.ans %in% c(1:10, 12:21, 22:25, 
            46:54)) {
            if (cont) {
                lower.aa <- min(yy/100)
                upper.aa <- max(yy * 100)
            }
            if (!cont) {
                lower.aa <- 0.1
                if (length(unique(fct3)) == 1) {
                  if (nth > 4) 
                    upper.aa <- exp(2 * sig.start)
                  else switch(nth, upper.aa <- exp(6 * sig.start), 
                    upper.aa <- exp(5 * sig.start), upper.aa <- exp(4 * 
                      sig.start), upper.aa <- exp(3 * sig.start))
                  if (quick.ans == 1) 
                    cat("\n\n Upper bound of parameter a was constrained based on\n                           the value", 
                      sig.start, "for sigma\n")
                }
                else {
                  lower.aa <- 1e-04
                  upper.aa <- 1000
                }
            }
            lower.bb <- -Inf
            upper.bb <- Inf
            if (model.ans %in% c(4:6, 9:10, 13:15, 16, 19:21, 
                23:25, 46:48, 50, 51:54)) 
                lower.bb <- 1e-06
            if (model.ans %in% 17:18) {
                if (increase == 1) {
                  lower.bb <- -Inf
                  upper.bb <- -1.01 * max(x)
                }
                else {
                  lower.bb <- eps
                  upper.bb <- Inf
                }
            }
            if (model.ans %in% 3) {
                if (increase == -1) {
                  lower.bb <- -Inf
                  upper.bb <- 0
                }
                if (increase == 1) {
                  lower.bb <- 0
                  upper.bb <- Inf
                }
            }
            if (model.ans %in% 18) {
                if (increase == 1) {
                  lower.bb <- -Inf
                  upper.bb <- 0
                }
                if (increase == -1) {
                  lower.bb <- 0
                  upper.bb <- Inf
                }
            }
            low.cc.default <- 1.3
            if (increase == 1) {
                lower.cc <- low.cc.default
                upper.cc <- 1e+06
            }
            if (increase == -1) {
                lower.cc <- 1e-06
                upper.cc <- 1/low.cc.default
            }
            if (increase == 0) {
                lower.cc <- 1e-06
                upper.cc <- 1000
            }
            if (incr.decr.no > 0) {
                lower.cc <- c(1e-18, 1e+18)
                upper.cc <- c(1e-18, 1e+18)
            }
            if (model.ans == 6 && ans.m6.sd == 2) {
                lower.cc <- 1e-06
                upper.cc <- 1000
            }
            if (model.ans == 49) {
                lower.cc <- 1e-06
                upper.cc <- 1000
            }
            if (incr.decr.no == 0) 
                if (cont && !MA.running && model.ans %in% c(4, 
                  5, 14, 15, 19, 20, 21, 24, 25, 46, 52, 54)) {
                  if (increase == 1) {
                    cc.feasible <- (1 + abs(CES)) * 1.01
                    if (cc.feasible > low.cc.default) {
                      lower.cc <- cc.feasible
                      cat("\nwarning: lower constraint on parameter cc has been increased to", 
                        lower.cc, "\ndue to value of CES\n")
                      cat("This might result in suboptimal fit of model", 
                        model.ans, "\n")
                    }
                  }
                  if (increase == -1) {
                    cc.feasible <- (1 - abs(CES))/1.01
                    if (upper.cc > cc.feasible) {
                      cat("\nwarning: upper constraint on parameter cc has been decreased to", 
                        cc.feasible, "\ndue to value of CES\n")
                      cat("This might result in suboptimal fit of model", 
                        model.ans, "\n")
                    }
                  }
                }
            constr.dd <- f.constr.dd(model.ans)
            lower.dd <- constr.dd[1]
            upper.dd <- constr.dd[2]
            if (model.ans == 1) {
                lower <- c(lower.aa)
                upper <- c(upper.aa)
            }
            if (model.ans %in% c(2, 7, 12, 17, 22)) {
                lower <- c(lower.aa, lower.bb)
                upper <- c(upper.aa, upper.bb)
            }
            if (model.ans %in% c(3, 8, 13, 18, 23, 51, 53)) {
                lower <- c(lower.aa, lower.bb, lower.dd)
                upper <- c(upper.aa, upper.bb, upper.dd)
            }
            if (model.ans %in% c(4, 9, 14, 19, 24)) {
                lower <- c(lower.aa, lower.bb, lower.cc)
                upper <- c(upper.aa, upper.bb, upper.cc)
            }
            if (model.ans %in% c(5, 6, 10, 15, 16, 20, 21, 25, 
                46, 49, 52, 54)) {
                lower <- c(lower.aa, lower.bb, lower.cc, lower.dd)
                upper <- c(upper.aa, upper.bb, upper.cc, upper.dd)
            }
            if (model.ans == 47) {
                lower.qq <- 5
                upper.qq <- 10
                lower <- c(lower.aa, lower.bb, lower.qq, lower.dd)
                upper <- c(upper.aa, upper.bb, upper.qq, upper.dd)
            }
            if (model.ans %in% c(48, 50)) {
                lower.aa.t <- lower.bb
                upper.aa.t <- upper.bb
                lower.bb.t <- 1e-06
                upper.bb.t <- 1e+06
                lower.dd.t <- 0.1
                upper.dd.t <- 10
                lower <- c(lower.aa, lower.cc, lower.dd, lower.aa.t, 
                  lower.bb.t, lower.dd.t)
                upper <- c(upper.aa, upper.cc, upper.dd, upper.aa.t, 
                  upper.bb.t, upper.dd.t)
            }
            if (cont && fit.ans == 1) {
                lower <- c(lower.var, lower)
                upper <- c(upper.var, upper)
            }
        }
        else if (model.ans == 11) {
            lower <- c(rep(lower.var, nr.var), regr.par)
            upper <- c(rep(upper.var, nr.var), regr.par)
        }
        else if (model.ans == 58) {
            lower.var <- lower[1]
            lower.FD.V <- lower[2]
            lower.ka <- lower[3]
            lower.ke <- lower[4]
            upper.var <- upper[1]
            upper.FD.V <- upper[2]
            upper.ka <- upper[3]
            upper.ke <- upper[4]
        }
        else if (model.ans == 59) {
            lower.var <- lower[1]
            lower.C1 <- lower[2]
            lower.C2 <- lower[3]
            lower.ka <- lower[4]
            lower.lab1 <- lower[5]
            lower.lab2 <- lower[6]
            upper.var <- upper[1]
            upper.C1 <- upper[2]
            upper.C2 <- upper[3]
            upper.ka <- upper[4]
            upper.lab1 <- upper[5]
            upper.lab2 <- upper[6]
        }
        else if (model.ans %in% c(31, 33, 57)) {
            lower.var <- lower[1]
            lower.aa <- lower[2]
            lower.b1 <- lower[3]
            lower.b2 <- lower[4]
            lower.c1 <- lower[5]
            lower.c2 <- lower[6]
            lower.dd <- lower[7]
            upper.var <- upper[1]
            upper.aa <- upper[2]
            upper.b1 <- upper[3]
            upper.b2 <- upper[4]
            upper.c1 <- upper[5]
            upper.c2 <- upper[6]
            upper.dd <- upper[7]
            if (model.ans == 33) {
                lower.d1 <- lower[7]
                lower.d2 <- lower[8]
                upper.d1 <- upper[7]
                upper.d2 <- upper[8]
            }
            if (model.ans == 57) {
                lower.d1 <- lower[6]
                lower.d2 <- lower[7]
                upper.d1 <- upper[6]
                upper.d2 <- upper[7]
            }
        }
        else if (model.ans == 65) {
            lower.cc <- lower[2]
            upper.cc <- upper[2]
            lower.k1 <- lower[3]
            upper.k1 <- upper[3]
            lower.ke <- lower[4]
            upper.ke <- upper[4]
        }
        else {
            if (cont) 
                dum <- 1
            else dum <- 0
            lower.var <- lower[dum]
            lower.aa <- lower[1 + dum]
            lower.bb <- lower[2 + dum]
            lower.cc <- lower[3 + dum]
            lower.dd <- lower[4 + dum]
            lower.ee <- lower[5 + dum]
            lower.ff <- lower[6 + dum]
            lower.gg <- lower[7 + dum]
            lower.hh <- lower[8 + dum]
            upper.var <- upper[dum]
            upper.aa <- upper[1 + dum]
            upper.bb <- upper[2 + dum]
            upper.cc <- upper[3 + dum]
            upper.dd <- upper[4 + dum]
            upper.ee <- upper[5 + dum]
            upper.ff <- upper[6 + dum]
            upper.gg <- upper[7 + dum]
            upper.hh <- upper[8 + dum]
            if (model.ans == 43) {
                lower.cc <- -1e+12
                upper.cc <- 1e+12
                lower[3 + dum] <- lower.cc
                upper[3 + dum] <- upper.cc
            }
        }
        if (tmp.quick) 
            if (length(xans) > 1) {
                lower <- c(lower, lower.RPF)
                upper <- c(upper, upper.RPF)
            }
        if (tmp.quick == F) {
            {
                vabcd <- f.text.par(ans.all, brief = T)
                if (incr.decr.no > 0) {
                }
                ans <- 1
                while (ans <= length(vabcd)) {
                  ans <- menu(c(paste(vabcd, ":     ", lower, 
                    " --- ", upper), "None, continue"), title = "\n For which parameter do you want to change the constraint? ")
                  if (ans <= length(vabcd)) {
                    cat("\n give infinite constraints as Inf or -Inf\n\n")
                    lower[ans] <- eval(parse(prompt = paste("give lower bound for ", 
                      vabcd[ans], "  ", "  > ")))
                    upper[ans] <- eval(parse(prompt = paste("give upper bound for ", 
                      vabcd[ans], "  > ")))
                  }
                }
                lower.var <- lower[1]
                lower.aa <- lower[2]
                lower.bb <- lower[3]
                upper.var <- upper[1]
                upper.aa <- upper[2]
                upper.bb <- upper[3]
                if (model.ans %in% c(4, 5, 6, 9, 10, 14, 15, 
                  19, 20, 21, 24, 25, 26, 30, 34, 43, 46, 49, 
                  52, 54, 60:63)) {
                  lower.cc <- lower[4]
                  upper.cc <- upper[4]
                  lower.dd <- lower[5]
                  upper.dd <- upper[5]
                  if (length(xans) > 1) {
                    lower.RPF <- lower[6]
                    upper.RPF <- upper[6]
                  }
                }
                if (model.ans %in% c(3, 8, 13, 18, 23, 51, 53)) {
                  lower.dd <- lower[4]
                  upper.dd <- upper[4]
                  if (length(xans) > 1) {
                    lower.RPF <- lower[5]
                    upper.RPF <- upper[5]
                  }
                }
                if (model.ans %in% c(2, 7, 12, 17, 22)) {
                  if (length(xans) > 1) {
                    lower.RPF <- lower[3]
                    upper.RPF <- upper[3]
                  }
                }
                if (model.ans %in% c(31, 33)) {
                  lower.b1 <- lower[3]
                  lower.b2 <- lower[4]
                  lower.c1 <- lower[5]
                  lower.c2 <- lower[6]
                  lower.dd <- lower[7]
                  upper.b1 <- upper[3]
                  upper.b2 <- upper[4]
                  upper.c1 <- upper[5]
                  upper.c2 <- upper[6]
                  upper.dd <- upper[7]
                  if (model.ans == 33) {
                    lower.d1 <- lower[7]
                    lower.d2 <- lower[8]
                    upper.d1 <- upper[7]
                    upper.d2 <- upper[8]
                  }
                  if (model.ans == 57) {
                    lower.d1 <- lower[6]
                    lower.d2 <- lower[7]
                    upper.d1 <- upper[6]
                    upper.d2 <- upper[7]
                  }
                }
                if (model.ans %in% c(48, 50)) {
                  lower.c <- lower[3]
                  lower.d <- lower[4]
                  lower.aa.t <- lower[5]
                  lower.bb.t <- lower[6]
                  lower.dd.t <- lower[7]
                  upper.c <- upper[3]
                  upper.d <- upper[4]
                  upper.aa.t <- upper[5]
                  upper.bb.t <- upper[6]
                  upper.dd.t <- upper[7]
                }
            }
        }
        if (model.ans == 1) {
            lb <- c(rep(lower.aa, nr.aa))
            ub <- c(rep(upper.aa, nr.aa))
        }
        if (model.ans %in% c(2, 7, 12, 17, 22, 35, 41, 44)) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb))
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb))
        }
        if (model.ans %in% c(3, 8, 13, 18, 23, 51, 53)) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
                rep(lower.dd, nr.dd))
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
                rep(upper.dd, nr.dd))
        }
        if (model.ans %in% c(4, 9, 14, 19, 24, 26, 27, 30, 55, 
            56)) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
                rep(lower.cc, nr.cc))
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
                rep(upper.cc, nr.cc))
        }
        if (model.ans %in% c(5, 6, 10, 15, 16, 20, 21, 25, 28, 
            29, 34, 37, 43, 49, 52, 54, 60:63)) {
            if (incr.decr.no > 0) 
                nr.cc <- 1
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
                rep(lower.cc, nr.cc), rep(lower.dd, nr.dd))
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
                rep(upper.cc, nr.cc), rep(upper.dd, nr.dd))
        }
        if (model.ans == 11) {
            lb <- regr.par
            ub <- regr.par
        }
        if (model.ans == 31) {
            lb <- c(rep(lower.aa, nr.aa), lower.b1, rep(lower.b2, 
                nr.bb), rep(lower.c1, nr.cc), lower.c2, lower.dd)
            ub <- c(rep(upper.aa, nr.aa), upper.b1, rep(upper.b2, 
                nr.bb), rep(upper.c1, nr.cc), upper.c2, upper.dd)
        }
        if (model.ans == 32) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
                lower.cc, lower.dd, lower.ee, lower.ff)
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
                upper.cc, upper.dd, upper.ee, upper.ff)
        }
        if (model.ans == 33) {
            lb <- c(rep(lower.aa, nr.aa), lower.b1, rep(lower.b2, 
                nr.bb), rep(lower.c1, nr.cc), lower.c2, lower.d1, 
                lower.d2)
            ub <- c(rep(upper.aa, nr.aa), upper.b1, rep(upper.b2, 
                nr.bb), rep(upper.c1, nr.cc), upper.c2, upper.d1, 
                upper.d2)
        }
        if (model.ans %in% 35:36) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb))
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb))
        }
        if (model.ans == 57) {
            nr.b1 <- nr.dd
            nr.b2 <- nr.bb
            lb <- c(rep(lower.aa, nr.aa), rep(lower.b1, nr.b1), 
                rep(lower.b2, nr.b2), rep(lower.c1, nr.cc), lower.d1, 
                lower.d2)
            ub <- c(rep(upper.aa, nr.aa), rep(upper.b1, nr.b1), 
                rep(upper.b2, nr.b2), rep(upper.c1, nr.cc), upper.d1, 
                upper.d2)
        }
        if (model.ans == 42) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, 2), lower.cc, 
                lower.dd)
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, 2), upper.cc, 
                upper.dd)
        }
        if (model.ans == 45) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
                lower.cc, lower.dd, lower.ee)
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
                upper.cc, upper.dd, upper.ee)
        }
        if (model.ans == 46) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
                lower.cc, lower.dd)
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
                upper.cc, upper.dd)
        }
        if (model.ans == 47) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
                lower.qq, lower.dd)
            ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
                upper.qq, upper.dd)
        }
        if (model.ans %in% c(48, 50)) {
            lb <- c(rep(lower.aa, nr.aa), rep(lower.cc, nr.cc), 
                rep(lower.dd, nr.dd), lower.aa.t, lower.bb.t, 
                lower.dd.t)
            ub <- c(rep(upper.aa, nr.aa), rep(upper.cc, nr.cc), 
                rep(upper.dd, nr.dd), upper.aa.t, upper.bb.t, 
                upper.dd.t)
        }
        if (model.ans == 58) {
            lb <- c(rep(lower.FD.V, nr.aa), rep(lower.ka, nr.bb), 
                rep(lower.ke, nr.cc))
            ub <- c(rep(upper.FD.V, nr.aa), rep(upper.ka, nr.bb), 
                rep(upper.ke, nr.cc))
        }
        if (model.ans == 59) {
            lb <- c(rep(lower.C1, nr.aa), rep(lower.C2, nr.bb), 
                lower.ka, lower.lab1, lower.lab2)
            ub <- c(rep(upper.C1, nr.aa), rep(upper.C2, nr.bb), 
                upper.ka, upper.lab1, upper.lab2)
        }
        if (model.ans == 64) {
            lb <- c(0, rep(0, nr.aa), rep(-Inf, nr.bb), rep(-Inf, 
                nr.cc), rep(-Inf, nr.dd))
            ub <- c(Inf, rep(Inf, nr.aa), rep(0, nr.bb), rep(0, 
                nr.cc), rep(0, nr.dd))
        }
        if (model.ans == 65) {
            lb <- c(rep(lower.cc, nr.aa), rep(lower.k1, nr.bb), 
                rep(lower.ke, nr.dd))
            ub <- c(rep(upper.cc, nr.aa), rep(upper.k1, nr.bb), 
                rep(upper.ke, nr.dd))
        }
        if (model.ans != 42) 
            if (!(model.ans == 11 && model.type == 1)) {
                if (length(xans) > 1) {
                  lb <- c(lb, rep(lower.RPF, nr.dosecol - 1))
                  ub <- c(ub, rep(upper.RPF, nr.dosecol - 1))
                }
            }
        if (cont && fit.ans == 1) {
            lb <- c(rep(lower.var, nr.var), lb)
            ub <- c(rep(upper.var, nr.var), ub)
        }
        if (!cont) {
            if (model.ans %in% c(2:3, 7:8)) 
                ub[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
            if (model.ans %in% c(4:5, 9:10)) 
                lb[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
            if (model.ans %in% c(17:21)) 
                lb[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
            if (model.ans %in% c(4, 5, 9, 10, 14, 15, 19, 20, 
                21, 24, 25, 46)) {
                lb[nr.aa + nr.bb + 1] <- 1e-06
                ub[nr.aa + nr.bb + 1] <- 0.999999
            }
            if (model.ans %in% c(12:15, 22:25, 51:54)) {
                if (nr.aa > 1 && nr.bb == 1) {
                  lb[2:(nr.aa + 1)] <- eps
                  ub[2:(nr.aa + 1)] <- Inf
                }
                else {
                  lb[(nr.aa + 1):(nr.aa + nr.bb)] <- eps
                  ub[(nr.aa + 1):(nr.aa + nr.bb)] <- Inf
                }
            }
            lb[nrp + 1] <- th.0.start
            ub[nrp + 1] <- th.0.start
            if (dtype == 6) 
                npar <- npar - max(fct3)
            if (nth > 1) {
                lb[(nrp + 2):(npar - 1)] <- -Inf
                ub[(nrp + 2):(npar - 1)] <- 1e-10
            }
            lb[npar] <- sig.start
            ub[npar] <- sig.start
            if (length(unique(fct3)) > 1) {
                lb[npar] <- sig.start/100
                ub[npar] <- sig.start * 100
            }
            if (0) {
                if (max(as.numeric(fct3)) > 1) {
                  lb[npar] <- 0
                  ub[npar] <- Inf
                }
            }
            if (dtype %in% c(4, 6)) 
                if (model.ans %in% c(3, 5, 18, 20, 13, 15, 23, 
                  25)) 
                  lb[npar - nth - 1] <- constr.steepness
            if (model.ans %in% c(4, 5, 9, 10, 14, 15, 24, 25)) {
                lb[nr.aa + nr.bb + 1] <- 0
                ub[nr.aa + nr.bb + 1] <- 1
            }
            if (dtype == 6) {
                lb <- c(rep(0, max(fct3)), lb)
                ub <- c(rep(Inf, max(fct3)), ub)
            }
        }
        ans.all$lb <- lb
        ans.all$ub <- ub
        ans.all$lower <- lower
        ans.all$upper <- upper
        text.par <- f.text.par(ans.all)
        if (cont || quick.ans > 1) 
            for (ii in (1:length(lb))) {
                if (lb[ii] == ub[ii]) 
                  text.par[ii] <- paste(text.par[ii], "(fixed)", 
                    sep = "")
            }
        ans.all$text.par <- text.par
        if (exists("track")) 
            print("f.constr.con:  END")
        return(ans.all)
    })
}



f.nlminb <- function(ans.all, tmp.quick = F) {
    if (exists("track2")) 
        print("f.nlminb")
    ans.all$tmp.quick <- tmp.quick
    with(ans.all, {
        th.par <- ans.all$th.par
        scale.dum <- abs(1/par.start)
        scale.dum[scale.dum == Inf] <- 1000
        if (dtype == 6) {
            nr.alfa <- max(fct3)
            scale.dum[1:nr.alfa] <- 1
        }
        else nr.alfa <- 0
        if (scale.ans && quick.ans == 1 && model.ans != 1 && 
            !((model.ans == 14 & model.type == 1) || (model.ans == 
                11 & model.type == 2))) {
            ans <- 1
            while (ans <= npar) {
                ans <- menu(c(paste(text.par, ":      ", signif(scale.dum, 
                  3)), "continue"), title = "\n For what parameter do you want to change scaling used in nlminb? ")
                if (ans <= npar) 
                  scale.dum[ans] <- eval(parse(prompt = paste("Give value larger than one ", 
                    text.par[ans], ": ")))
            }
        }
        ans.all$scale.dum <- scale.dum
        trace.tmp <- F
        if (cont) {
            if (incr.decr.no == 0) 
                if (model.ans %in% c(4, 5, 9, 10, 14, 15, 19, 
                  20, 24, 25)) {
                  cc.start <- par.start[(nr.var + nr.aa + nr.bb + 
                    1):(nr.var + nr.aa + nr.bb + nr.cc)]
                  if (increase == -1 && cc.start >= (1 - abs(CES))) 
                    cc.start <- 0.9 * (1 - abs(CES))
                  if (increase == 1 && cc.start <= (1 + abs(CES))) 
                    cc.start <- 1.1 * (1 + abs(CES))
                  par.start[(nr.var + nr.aa + nr.bb + 1):(nr.var + 
                    nr.aa + nr.bb + nr.cc)] <- cc.start
                }
            loglik.check <- f.lik.con(par.start, x, y, dtype, 
                fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                ttt = ttt, trace.tmp = trace.tmp, fct4 = fct4, 
                fct5 = fct5, cens.up = cens.up, par.tmp = NA, 
                increase = increase, x.mn = x.mn, ref.lev = ref.lev, 
                sign.q = sign.q, ans.m6.sd = ans.m6.sd, Mx = Mx, 
                x1 = x1, x2 = x2, cc.inf = cc.inf, incr.decr.no = incr.decr.no)
            count <- 0
            while (!is.finite(loglik.check) && model.ans %in% 
                c(17, 18, 22, 23) && count < 30) {
                ced.start <- par.start[(nr.var + nr.aa + 1):(nr.var + 
                  nr.aa + nr.bb)]
                ced.start <- 1.5 * ced.start
                par.start[(nr.var + nr.aa + 1):(nr.var + nr.aa + 
                  nr.bb)] <- ced.start
                loglik.check <- f.lik.con(par.start, x, y, dtype, 
                  fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                  nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                  ttt = ttt, trace.tmp = trace.tmp, fct4 = fct4, 
                  fct5 = fct5, cens.up = cens.up, par.tmp = NA, 
                  increase = increase, x.mn = x.mn, ref.lev = ref.lev, 
                  sign.q = sign.q, ans.m6.sd = ans.m6.sd, Mx = Mx, 
                  x1 = x1, x2 = x2, cc.inf = cc.inf, incr.decr.no = incr.decr.no)
                count <- count + 1
            }
            if (model.ans == -51) {
                cat("\n")
                print("f.nlminb start")
                print(model.ans)
                print(CES)
                print(increase)
                print(cbind(text.par, par.start, lb, ub, scale.dum))
                print(loglik.check)
                f.press.key.to.continue()
            }
            fit.res <- nlminb(par.start, f.lik.con, scale = scale.dum, 
                lower = lb, upper = ub, control = lst.control, 
                x = x, y = y, dtype = dtype, fct1 = fct1, fct2 = fct2, 
                fct3 = fct3, fct4 = fct4, fct5 = fct5, model.ans = model.ans, 
                mn.log = mn.log, sd2.log = sd2.log, nn = nn, 
                Vdetlim = Vdetlim, CES = CES, twice = twice, 
                ttt = ttt, trace.tmp = F, cens.up = cens.up, 
                lb = lb, ub = ub, par.tmp = par.start, increase = increase, 
                x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                ans.m6.sd = ans.m6.sd, Mx = Mx, x1 = x1, x2 = x2, 
                cc.inf = cc.inf, incr.decr.no = incr.decr.no)
            if (is.na(fit.res$obj)) 
                fit.res$obj <- 1e-12
            count <- 0
            while (is.na(fit.res$par[1]) || (fit.res$obj == 0) || 
                abs(fit.res$obj) == 1e-12 || !is.finite(fit.res$obj) && 
                count <= 20) {
                cat("\n      model is refitted with other scaling parameter")
                scale.dum <- rep(0.5, length(par.start))
                scale.dum <- 2 * scale.dum
                fit.res <- nlminb(par.start, f.lik.con, scale = scale.dum, 
                  lower = lb, upper = ub, control = lst.control, 
                  x = x, y = y, dtype = dtype, fct1 = fct1, fct2 = fct2, 
                  fct3 = fct3, fct4 = fct4, fct5 = fct5, model.ans = model.ans, 
                  mn.log = mn.log, sd2.log = sd2.log, nn = nn, 
                  Vdetlim = Vdetlim, CES = CES, twice = T, cens.up = cens.up, 
                  lb = lb, ub = ub, par.tmp = fit.res$par, increase = increase, 
                  x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                  ans.m6.sd = ans.m6.sd, Mx = Mx, x1 = x1, x2 = x2, 
                  cc.inf = cc.inf, incr.decr.no = incr.decr.no)
                if (is.na(fit.res$obj)) 
                  fit.res$obj <- 1e-12
                count <- count + 1
            }
            ans.all$MLE <- fit.res$par
            if (model.ans == 47 || (model.ans == 6 && ans.m6.sd == 
                2)) 
                ans.all$regr.par <- ans.all$MLE
            else ans.all$regr.par <- ans.all$MLE[-(1:max(fct3))]
            if (model.ans %in% c(31, 33, 57)) 
                ans.all$regr.par.matr <- f.pars.m31(ans.all)$regr.par.matr
            else ans.all$regr.par.matr <- f.pars(ans.all)$regr.par.matr
            if (0) {
                loglik.tmp <- -f.lik.con(MLE, x, y, dtype, fct1, 
                  fct2, fct3, model.ans, mn.log, sd2.log, nn, 
                  Vdetlim = Vdetlim, CES = CES, twice = twice, 
                  ttt = ttt, trace.tmp = T, fct4 = fct4, fct5 = fct5, 
                  cens.up = cens.up, par.tmp = NA, increase = increase, 
                  x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                  ans.m6.sd = ans.m6.sd, Mx = Mx, x1 = x1, x2 = x2, 
                  cc.inf = cc.inf, incr.decr.no = incr.decr.no)
                print(loglik.tmp)
                f.press.key.to.continue()
            }
        }
        if (!cont) {
            loglik.check <- f.lik.cat(par.start, x = x, y = y, 
                kk = kk, nn = nn, dtype = dtype, fct1 = fct1, 
                fct2 = fct2, nrp = nrp, nth = nth, th.par = th.par, 
                nr.aa = nr.aa, nr.bb = nr.bb, model.ans = model.ans, 
                model.type = model.type, ans.nobg = 0, CES = CES, 
                CES.cat = CES.cat, ttt = ttt, twice = twice, 
                cens.up = cens.up, fct3 = fct3, fct4 = fct4, 
                fct5 = fct5, x.full = x.full, fct1.full = fct1.full, 
                fct2.full = fct2.full, trace.tmp = trace.tmp, 
                ces.ans = ces.ans, cc.inf = cc.inf, decr.zz = decr.zz, 
                CES1 = CES1, CES2 = CES2, nn.tot = nn.tot, kk.tot = kk.tot, 
                ref.lev = ref.lev, x1 = x1, x2 = x2, output = output, 
                Mx = Mx)
            if (model.ans == -5) {
                cat("\n")
                print("f.nlminb start")
                print(model.ans)
                print(CES)
                print(CES.cat)
                print(cbind(text.par, par.start, lb, ub, scale.dum))
                print(loglik.check)
                f.press.key.to.continue()
            }
            retry <- FALSE
            if (!is.finite(loglik.check) || abs(loglik.check) == 
                2e+10) 
                retry <- TRUE
            if (retry) {
                count <- 0
                while (retry && (model.type == 2 && model.ans %in% 
                  c(13, 18, 23)) && count < 15) {
                  ced.start <- par.start[(nr.alfa + nr.aa + 1):(nr.alfa + 
                    nr.aa + nr.bb)]
                  ced.start <- 1.5 * ced.start
                  par.start[(nr.alfa + nr.aa + 1):(nr.alfa + 
                    nr.aa + nr.bb)] <- ced.start
                  if (count < 3) 
                    cat("\nModel re-fitted with new start values\n")
                  loglik.check <- f.lik.cat(par.start, x = x, 
                    y = y, kk = kk, nn = nn, dtype = dtype, fct1 = fct1, 
                    fct2 = fct2, nrp = nrp, nth = nth, th.par = th.par, 
                    nr.aa = nr.aa, nr.bb = nr.bb, model.ans = model.ans, 
                    model.type = model.type, ans.nobg = 0, CES = CES, 
                    CES.cat = CES.cat, ttt = ttt, twice = twice, 
                    cens.up = cens.up, fct3 = fct3, fct4 = fct4, 
                    fct5 = fct5, x.full = x.full, fct1.full = fct1.full, 
                    fct2.full = fct2.full, trace.tmp = trace.tmp, 
                    ces.ans = ces.ans, cc.inf = cc.inf, decr.zz = decr.zz, 
                    CES1 = CES1, CES2 = CES2, nn.tot = nn.tot, 
                    kk.tot = kk.tot, ref.lev = ref.lev, Mx = Mx, 
                    x1 = x1, x2 = x2, output = output)
                  if (is.finite(loglik.check) && abs(loglik.check) < 
                    2e+10) 
                    retry <- FALSE
                  count <- count + 1
                }
            }
            fit.res <- nlminb(par.start, f.lik.cat, scale = scale.dum, 
                lower = lb, upper = ub, control = lst.control, 
                x = x, y = y, kk = kk, nn = nn, dtype = dtype, 
                fct1 = fct1, fct2 = fct2, nrp = nrp, nth = nth, 
                th.par = th.par, nr.aa = nr.aa, nr.bb = nr.bb, 
                model.ans = model.ans, model.type = model.type, 
                ans.nobg = 0, CES = CES, CES.cat = CES.cat, ttt = ttt, 
                twice = twice, cens.up = cens.up, fct3 = fct3, 
                fct4 = fct4, fct5 = fct5, x.full = x.full, fct1.full = fct1.full, 
                fct2.full = fct2.full, trace.tmp = trace.tmp, 
                ces.ans = ces.ans, cc.inf = cc.inf, decr.zz = decr.zz, 
                CES1 = CES1, CES2 = CES2, nn.tot = nn.tot, kk.tot = kk.tot, 
                ref.lev = ref.lev, Mx = Mx, x1 = x1, x2 = x2, 
                output = output)
            if (dtype == 6 && model.type == 2 && model.ans %in% 
                c(13, 15, 23, 25) && quick.ans > 1) {
                par.start <- f.start.cat(ans.all, tmp.quick = TRUE)$par.start
                scale.dum <- abs(1/par.start)
                scale.dum[scale.dum == Inf] <- 1000
                scale.dum[1:nr.alfa] <- 1
                fit.res.test <- nlminb(par.start, f.lik.cat, 
                  scale = scale.dum, lower = lb, upper = ub, 
                  control = lst.control, x = x, y = y, kk = kk, 
                  nn = nn, dtype = dtype, fct1 = fct1, fct2 = fct2, 
                  nrp = nrp, nth = nth, th.par = th.par, nr.aa = nr.aa, 
                  nr.bb = nr.bb, model.ans = model.ans, model.type = model.type, 
                  ans.nobg = 0, CES = CES, CES.cat = CES.cat, 
                  ttt = ttt, twice = twice, cens.up = cens.up, 
                  fct3 = fct3, fct4 = fct4, fct5 = fct5, x.full = x.full, 
                  fct1.full = fct1.full, fct2.full = fct2.full, 
                  trace.tmp = trace.tmp, ces.ans = ces.ans, cc.inf = cc.inf, 
                  decr.zz = decr.zz, CES1 = CES1, CES2 = CES2, 
                  nn.tot = nn.tot, kk.tot = kk.tot, ref.lev = ref.lev, 
                  Mx = Mx, x1 = x1, x2 = x2, output = output)
                if (fit.res.test$objective < fit.res$objective) 
                  fit.res <- fit.res.test
            }
            if (dtype == 3) 
                if (model.type == 2 && model.ans %in% c(12:16, 
                  22:25, 46:47)) {
                  ced.start <- par.start[(nr.aa + 1):(nr.aa + 
                    nr.bb)]
                  par.start[(nr.aa + 1):(nr.aa + nr.bb)] <- ced.start/10
                  fit.res.test <- nlminb(par.start, f.lik.cat, 
                    scale = scale.dum, lower = lb, upper = ub, 
                    control = lst.control, x = x, y = y, kk = kk, 
                    nn = nn, dtype = dtype, fct1 = fct1, fct2 = fct2, 
                    nrp = nrp, nth = nth, th.par = th.par, nr.aa = nr.aa, 
                    nr.bb = nr.bb, model.ans = model.ans, model.type = model.type, 
                    ans.nobg = 0, CES = CES, CES.cat = CES.cat, 
                    ttt = ttt, twice = twice, cens.up = cens.up, 
                    fct3 = fct3, fct4 = fct4, fct5 = fct5, x.full = x.full, 
                    fct1.full = fct1.full, fct2.full = fct2.full, 
                    trace.tmp = trace.tmp, ces.ans = ces.ans, 
                    cc.inf = cc.inf, decr.zz = decr.zz, CES1 = CES1, 
                    CES2 = CES2, nn.tot = nn.tot, kk.tot = kk.tot, 
                    ref.lev = ref.lev, Mx = Mx, x1 = x1, x2 = x2, 
                    output = output)
                  if (fit.res.test$objective < fit.res$objective) 
                    fit.res <- fit.res.test
                }
            if (model.ans == -10) {
                print(" f.nlminb cont")
                print(model.ans)
                print(CES)
                print(CES.cat)
                print(cbind(text.par, par.start, lb, ub, scale.dum))
                print(fit.res$objective)
                print(fit.res$par)
                f.press.key.to.continue()
            }
            count <- 0
            while (is.na(fit.res$par[1]) || (fit.res$obj == 1e-12) | 
                is.na(fit.res$obj)) {
                cat("\n     model", model.ans, "refitted with other scaling parameter")
                if (model.ans == 14 && model.type == 1 && dtype == 
                  6) 
                  par.start[1] <- eval(parse(prompt = paste("give start value for alfa", 
                    "  > ")))
                scale.dum <- rep(0.5, length(par.start))
                scale.dum <- 2 * scale.dum
                count <- count + 1
                fit.res <- nlminb(par.start, f.lik.cat, scale = scale.dum, 
                  lower = lb, upper = ub, control = lst.control, 
                  x = x, y = y, kk = kk, nn = nn, dtype = dtype, 
                  fct1 = fct1, fct2 = fct2, nrp = nrp, nth = nth, 
                  th.par = th.par, nr.aa = nr.aa, nr.bb = nr.bb, 
                  model.ans = model.ans, model.type = model.type, 
                  ans.nobg = 0, CES = CES, CES.cat = CES.cat, 
                  ttt = ttt, twice = twice, cens.up = cens.up, 
                  fct3 = fct3, fct5 = fct5, x.full = x.full, 
                  fct1.full = fct1.full, fct2.full = fct2.full, 
                  cc.inf = cc.inf, trace.tmp = trace.tmp, ces.ans = ces.ans, 
                  decr.zz = decr.zz, CES1 = CES1, CES2 = CES2, 
                  ref.lev = ref.lev, Mx = Mx, x1 = x1, x2 = x2, 
                  output = output)
                if (count > 50) {
                  fit.res$obj <- -1e+10
                  fit.res$par <- 0
                }
            }
            ans.all$MLE <- fit.res$par
            if (model.type == 1) 
                ans.all$regr.par <- ans.all$MLE
            if (model.type == 2) 
                if (length(unique(fct3)) == 1) {
                  par.lst <- f.split.par(ans.all$MLE, nrp, nth, 
                    dtype, fct3)
                  ans.all$regr.par <- par.lst$regr.par
                  ans.all$th.par <- par.lst$th.par
                  ans.all$sig.par <- par.lst$sig.par
                  ans.all <- f.pars(ans.all)
                }
                else {
                  ans.all$get.gr.txt <- TRUE
                  ans.all$regr.par <- ans.all$MLE[1:nrp]
                  ans.all <- f.pars.frq(ans.all)
                }
        }
        ans.all$loglik <- round(-fit.res$objective, 2)
        if (!tmp.quick) 
            if (quick.ans == 1 && length(ans.all$loglik.start) > 
                0) 
                if (abs(ans.all$loglik - loglik.start) < 0.1) {
                  cat("\n\nATTENTION:  log-likelihood hardly changed\n\n")
                  f.press.key.to.continue()
                }
        if ((model.type == 1 && model.ans == 14) || (cont && 
            model.ans == 11)) 
            npar.aic <- npar
        else npar.aic <- length(ans.all$MLE) - sum(lb == ub)
        ans.all$aic <- 2 * npar.aic - 2 * ans.all$loglik
        ans.all$npar.aic <- npar.aic
        if (!is.finite(ans.all$loglik)) 
            ans.all$loglik <- -1e+12
        mess <- fit.res$message
        ans.all$converged <- f.converged(mess, fit.res$conv, 
            tmp.quick = tmp.quick)
        if (dtype == 6 && (model.ans == 14 & model.type == 1)) {
            dum <- list()
            dum$loglik <- ans.all$loglik
            dum$MLE <- ans.all$MLE
            ans.all$full.model <- dum
            alfa.start <- ans.all$MLE[1]
            nn.dum <- table(x, x)
            nn.dum <- nn.dum[nn.dum != 0]
            if (length(nn.dum) > 15 && alfa.start > 20) 
                alfa.start <- 5
            ans.all$alfa.start <- alfa.start
        }
        ans.all$fitted <- T
        ans.all$l.ty <- 1
        ans.all$new.model <- F
        if (0) {
            ans.all$fit.res <- fit.res
            f.vcov(ans.all)
        }
        if (exists("track2")) 
            print("f.nlminb : END")
        return(ans.all)
    })
}



f.pars <- function(ans.all) {
    if (exists("track2")) 
        print("f.pars")
    if (ans.all$fit.ans == 1) 
        if (ans.all$model.ans == 47 || (ans.all$model.ans == 
            6 && ans.all$ans.m6.sd == 2)) 
            return(f.pars.m6(ans.all))
    with(ans.all, {
        if (!cont && model.ans == 14 && model.type == 1) {
            ans.all$regr.par.matr <- 0
            return(ans.all)
        }
        else if (cont && model.ans == 11) {
            ans.all$regr.par.matr <- 0
            return(ans.all)
        }
        gr.txt.save <- ans.all$gr.txt
        nrp <- length(regr.par)
        regr.par.matr <- numeric()
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        nr.cc <- max(fct4)
        nr.dd <- max(fct5)
        nr.var <- max(fct3)
        if (dtype == 4 && max(fct3) > 1) 
            nr.var <- 1
        if (length(ans.all$nr.cc) == 0) 
            nr.cc <- 1
        if (length(ans.all$nr.dd) == 0) 
            nr.dd <- 1
        if (length(ans.all$nr.var) == 0) 
            nr.var <- 1
        if (model.ans == 36) {
            nr.cc <- 0
            nr.dd <- 0
        }
        if (nr.aa == 1 && nr.bb == 1 && nr.cc == 1 && nr.dd == 
            1 && nr.var < 2) {
            ans.all$regr.par.matr <- matrix(regr.par, nrow = 1)
            ans.all$nr.gr <- 1
            ans.all$gr.txt <- ""
            if (max(fct3) > 1) {
                ans.all$gr.txt <- fct3.txt
                ans.all$regr.par.matr <- matrix(regr.par, byrow = TRUE, 
                  nrow = length(fct3.txt), ncol = length(regr.par))
            }
            if (length(xans) > 1) 
                ans.all$gr.txt <- gr.txt.save
            if (exists("track2")) 
                print("end of f.pars")
            return(ans.all)
        }
        if (cont || model.type == 2) {
            if (model.ans %in% c(3, 8, 13, 18, 23, 51, 53)) 
                nr.cc <- 0
            if (model.ans %in% c(4, 9, 14, 19, 24)) 
                nr.dd <- 0
            if (model.ans %in% c(2, 12, 17, 22)) {
                nr.cc <- 0
                nr.dd <- 0
            }
            cc <- regr.par[nr.aa + nr.bb + 1]
            dd <- regr.par[nr.aa + nr.bb + nr.cc + 1]
            if (cont && model.ans == 31) 
                dd <- regr.par[(nr.aa + nr.bb + 1):(nr.aa + nr.bb + 
                  4)]
            if (cont && model.ans == 32) 
                dd <- regr.par[(nr.aa + nr.bb + 1):(nr.aa + nr.bb + 
                  4)]
            if (cont && model.ans == 33) 
                dd <- regr.par[(nr.aa + nr.bb + 1):(nr.aa + nr.bb + 
                  5)]
            if (cont && model.ans == 45) 
                dd <- regr.par[(nr.aa + nr.bb + 1):(nr.aa + nr.bb + 
                  3)]
        }
        if (cont && model.ans %in% c(35, 41, 44)) 
            nr.cc <- 0
        if (cont && model.ans %in% c(26, 27, 30, 35, 41, 44, 
            55, 56)) 
            nr.dd <- 0
        if (!cont && model.type == 1) {
            if (model.ans %in% c(7, 20)) {
                nr.dd <- 1
                dd <- regr.par[nr.aa + nr.bb + nr.cc + 1]
            }
            else nr.dd <- 0
        }
        nr.subgr <- max(nr.aa, nr.bb, nr.cc, nr.dd)
        if (nr.aa == 1) 
            fct1.txt <- rep("", nr.subgr)
        if (nr.bb == 1) 
            fct2.txt <- rep("", nr.subgr)
        if (nr.cc == 1) 
            fct4.txt <- rep("", nr.subgr)
        if (nr.dd == 1) 
            fct5.txt <- rep("", nr.subgr)
        if (identical(fct1, fct2)) 
            fct2.txt <- rep("", nr.subgr)
        if (identical(fct1, fct4)) 
            fct4.txt <- rep("", nr.subgr)
        if (identical(fct4, fct2)) 
            fct4.txt <- rep("", nr.subgr)
        gr.txt <- character(0)
        if (!cont && model.type == 2 && (model.ans %in% c(12:15, 
            22:25)) && ces.ans %in% c(1:3, 5)) {
            if (nrp > (nr.aa + nr.bb)) 
                par.rest <- regr.par[(nr.aa + nr.bb + 1):nrp]
            else par.rest <- numeric()
            if (nr.aa > 1 && nr.bb == 1) 
                for (ii in 1:nr.aa) {
                  par.tmp <- c(regr.par[1], regr.par[ii + 1], 
                    par.rest)
                  regr.par.matr <- rbind(regr.par.matr, par.tmp)
                }
            if (nr.aa == 1 & nr.bb > 1) 
                for (jj in 1:nr.bb) {
                  par.tmp <- c(regr.par[1], regr.par[jj + 1], 
                    par.rest)
                  regr.par.matr <- rbind(regr.par.matr, par.tmp)
                }
            if (nr.aa > 1 & nr.bb > 1) 
                for (jj in 1:nr.bb) {
                  par.tmp <- c(regr.par[jj], regr.par[jj + nr.aa], 
                    par.rest)
                  regr.par.matr <- rbind(regr.par.matr, par.tmp)
                }
            ans.all$regr.par.matr <- regr.par.matr
            if (quick.ans > 1 && length(ans.all$covar.txt) > 
                0) 
                fct1.txt <- covar.txt
            kk <- 1
            for (jj in 1:nr.bb) {
                f1 <- fct1[fct2 == jj]
                f1.lev <- levels(factor(f1))
                for (ii.index in f1.lev) {
                  ii <- as.numeric(ii.index)
                  if (!twice) 
                    gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                      sep = "-")
                  if (twice) {
                    gr.txt[kk] <- paste(fct1.txt[kk], sep = "-")
                  }
                  kk <- kk + 1
                }
            }
            ans.all$nr.gr <- length(ans.all$regr.par.matr[, 1])
            ans.all$gr.txt <- gr.txt
            if (length(xans) > 1) 
                ans.all$gr.txt <- gr.txt.save
            if (exists("track2")) 
                print("f.pars LVM: END ")
            return(ans.all)
        }
        if (nr.dd < 2) {
            if (length(fct4) == 1) 
                fct4 <- rep(1, length(x))
            par.tmp <- rep(NA, length(regr.par))
            kk <- 1
            for (jj in 1:nr.bb) {
                f1 <- fct1[fct2 == jj]
                f1.lev <- levels(factor(f1))
                for (ii.index in f1.lev) {
                  ii <- as.numeric(ii.index)
                  f4 <- fct4[fct1 == ii & fct2 == jj]
                  f4.lev <- levels(factor(f4))
                  for (mm.index in f4.lev) {
                    mm <- as.numeric(mm.index)
                    if (cont && model.ans %in% c(31, 32, 33, 
                      45)) {
                      mm <- 0
                      nr.cc <- 0
                    }
                    par.tmp <- c(regr.par[ii], regr.par[nr.aa + 
                      jj])
                    if (nr.cc > 0) 
                      par.tmp <- c(par.tmp, regr.par[nr.aa + 
                        nr.bb + mm])
                    if (nr.dd > 0) 
                      par.tmp <- c(par.tmp, dd)
                    if (length(xans) > 1) {
                      RPF.vec <- ans.all$regr.par[(length(ans.all$regr.par) - 
                        nr.dosecol + 2):length(ans.all$regr.par)]
                      par.tmp <- c(par.tmp, RPF.vec)
                    }
                    regr.par.matr <- rbind(regr.par.matr, par.tmp)
                    gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                      fct4.txt[mm], sep = "-")
                    kk <- kk + 1
                  }
                }
            }
            n.col <- length(par.tmp)
            regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                ncol = n.col)
        }
        if (exists("track2")) 
            print("f.pars, continue 2")
        if (nr.cc < 2 && nr.dd > 1) {
            if (length(fct5) == 1) 
                fct5 <- rep(1, length(x))
            par.tmp <- rep(NA, length(regr.par))
            kk <- 1
            for (jj in 1:nr.bb) {
                f1 <- fct1[fct2 == jj]
                f1.lev <- levels(factor(f1))
                for (ii.index in f1.lev) {
                  ii <- as.numeric(ii.index)
                  f5 <- fct5[fct1 == ii & fct2 == jj]
                  f5.lev <- levels(factor(f5))
                  for (mm.index in f5.lev) {
                    mm <- as.numeric(mm.index)
                    par.tmp <- c(regr.par[ii], regr.par[nr.aa + 
                      jj])
                    if (nr.cc > 0) 
                      par.tmp <- c(par.tmp, regr.par[nr.aa + 
                        nr.bb + 1])
                    if (nr.dd > 0) 
                      par.tmp <- c(par.tmp, regr.par[nr.aa + 
                        nr.bb + nr.cc + mm])
                    if (length(xans) > 1) {
                      RPF.vec <- ans.all$regr.par[(length(ans.all$regr.par) - 
                        nr.dosecol + 2):length(ans.all$regr.par)]
                      par.tmp <- c(par.tmp, RPF.vec)
                    }
                    regr.par.matr <- rbind(regr.par.matr, par.tmp)
                    gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                      fct5.txt[mm], sep = "-")
                    kk <- kk + 1
                  }
                }
            }
            n.col <- length(par.tmp)
            regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                ncol = n.col)
        }
        if (exists("track2")) 
            print("f.pars, continue 3")
        ans.all$regr.par.matr <- regr.par.matr
        ans.all$gr.txt <- gr.txt
        ans.all$nr.gr <- length(regr.par.matr[, 1])
        if (length(xans) > 1) 
            ans.all$gr.txt <- gr.txt.save
        if (exists("track2")) 
            print("f.pars: END ")
        return(ans.all)
    })
}



f.converged <- function(mess, conv.out, tmp.quick = F) {
    if (exists("track2")) 
        print("f.converged")
    converged <- F
    if (mode(mess) == "NULL") {
        converged <- F
        if (!tmp.quick) 
            cat("\nATTENTION: convergence message is NULL\n")
    }
    else converged <- 1 - conv.out
    if (!tmp.quick) 
        cat("\n CONVERGION TYPE: ", mess, "\n")
    if (exists("track2")) 
        print("f.converged:  END")
    return(converged)
}



f.bb.con <- function(model.ans, cc, dd, CED, CES, ref.lev = NA, CED.ref = NA, 
    incr.decr.no = 0, cont, aa = NA) {
    if (exists("track2")) 
        print("f.bb.con")
    if (model.ans == 12) 
        bb <- logb(CES + 1)/CED
    if (model.ans == 14) 
        bb <- -logb((CES + 1 - cc)/(1 - cc))/CED
    if (model.ans %in% c(13, 15, 50)) {
        if (incr.decr.no > 0) 
            CES.tmp <- CES
        else {
            if ((cc[1] < 1)) 
                CES.tmp <- -abs(CES)
            else CES.tmp <- abs(CES)
        }
        bb <- CED/(-log(1 - log(CES.tmp + 1)/log(cc)))^(1/dd)
    }
    if (model.ans == 22) 
        bb <- -CED * (1 + CES)/CES
    if (model.ans %in% c(23, 25)) {
        if (incr.decr.no > 0) 
            CES.tmp <- CES
        else {
            if ((cc[1] < 1)) 
                CES.tmp <- -abs(CES)
            else CES.tmp <- abs(CES)
        }
        dum <- log(CES.tmp + 1)/log(cc)
        bb <- CED * ((1 - dum)/dum)^(1/dd)
    }
    if (model.ans == 24) 
        bb <- CED/(CES/(cc - 1 - CES))
    if (model.ans == 46) {
        if ((cc[1] < 1)) 
            CES.tmp <- -abs(CES)
        else CES.tmp <- abs(CES)
        CED <- CED.ref/CED
        CED[ref.lev] <- CED.ref
        if (cont) 
            bb <- CED/(-log(1 - log(CES.tmp + 1)/log(cc)))^(1/dd)
    }
    if (model.ans %in% c(51, 52)) {
        if (incr.decr.no > 0) 
            CES.tmp <- CES
        else {
            if ((cc[1] < 1)) 
                CES.tmp <- -abs(CES)
            else CES.tmp <- abs(CES)
        }
        dum <- log(CES.tmp + 1)/log(cc)
        dum <- (-log(dum))^(1/dd)
        bb <- CED * dum
    }
    if (model.ans %in% c(53, 54)) {
        if (incr.decr.no > 0) 
            CES.tmp <- CES
        else {
            if ((cc[1] < 1)) 
                CES.tmp <- -abs(CES)
            else CES.tmp <- abs(CES)
        }
        dum <- log(CES.tmp + 1)/log(cc)
        log.bb <- qnorm(dum) - dd * log(CED)
        bb <- exp(log.bb)
    }
    if (model.ans == 60) {
        dum <- log(1 - log(CES + 1)/(aa * (cc - 1)))
        bb <- -dum/(CED^dd)
    }
    if (model.ans == 61) {
        dum <- log(1 - log(CES + 1)/(aa * (cc - 1)))
        bb <- -dum/(CED^dd)
    }
    if (model.ans == 63) {
        dum <- qnorm((1/aa) * log(CES + 1) + pnorm(cc))
        bb <- (dum - cc)/(CED^dd)
    }
    if (exists("track2")) 
        print("f.bb.con:END")
    return(bb)
}



f.grubb <- function(ss = 25, alfa = 0.05) {
    alfa <- alfa/2
    t.crit <- qt(alfa/ss, ss - 2)
    G <- (ss - 1)/sqrt(ss)
    G <- G * sqrt(t.crit^2/(ss - 2 + t.crit^2))
    return(G)
}



f.CI <- function(ans.all, display_plots = TRUE) {
    if (exists("track")) 
        print("f.CI")
    with(ans.all, {
        if (quick.ans == 1) 
            ans.all$trace.plt <- T
        profile.out <- f.profile.all(ans.all, display_plots = display_plots)
        conf.int <- profile.out$conf.int
        update <- F
        count <- 0
        loglik.tmp <- profile.out$loglik - 0.5
        ans.all$CED.low <- profile.out$CED.low
        ans.all$loglik.low <- profile.out$loglik.low
        ans.all$CED.upp <- profile.out$CED.upp
        ans.all$loglik.upp <- profile.out$loglik.upp
        while (sum(profile.out$MLE.new) != 0) {
            count <- count + 1
            if (interrupt) {
                cat(" \n\n local optimum found .... calculation of CI is re-started\n")
                cat(" improved fit will be stored in R-object called newfit\n\n")
            }
            if (profile.out$loglik < loglik.tmp - 0.02) {
                cat("\n\nATTENTION: log-likelihood fluctuates, CI is not calculated\n\n")
                profile.out$conf.int <- matrix(NA, ncol = 2)
                profile.out$MLE.new <- 0
            }
            else {
                ans.all$par.start <- profile.out$MLE.new
                ans.all <- f.nlminb(ans.all)
                cat("\n\n new log-likelihood: ", ans.all$loglik, 
                  "\n\n")
                if (!WAPP && !gui) 
                  f.store.results(ans.all, "newfit")
                profile.out <- f.profile.all(ans.all, display_plots = display_plots)
                loglik.tmp <- profile.out$loglik
                update <- T
                if (count > 50) {
                  cat("\n\nno global optimum found\n\n")
                  profile.out$MLE.new <- 0
                  profile.out$conf.int <- matrix(NA, ncol = 2)
                }
                conf.int <- profile.out$conf.int
            }
        }
        if (update) {
            if (cont) 
                ans.all$regr.par <- ans.all$MLE[-(1:max(fct3))]
            if (!cont) {
                if (dtype == 6) 
                  nr.alfa <- max(fct3)
                else nr.alfa <- 0
                CED <- ans.all$MLE[(nr.alfa + nr.var + nr.aa + 
                  1):(nr.alfa + nr.var + nr.aa + nr.bb)]
                CED.matr <- matrix(CED, ncol = 1)
                ans.all$CED <- signif(CED, 2)
                ans.all$CED.matr <- CED.matr
                if (model.type == 1) {
                  ans.all$regr.par <- ans.all$MLE[1:nrp]
                  if (dtype == 6) 
                    ans.all$regr.par <- ans.all$regr.par[-1]
                }
                if (model.type == 2) {
                  par.lst <- f.split.par(ans.all$MLE, nrp, nth, 
                    dtype, fct3)
                  ans.all$regr.par <- par.lst$regr.par
                  ans.all$th.par <- par.lst$th.par
                  ans.all$sig.par <- par.lst$sig.par
                }
                ans.all$show <- f.show.cat(ans.all)
            }
        }
        low <- conf.int[, 1]
        upp <- conf.int[, 2]
        lst1 <- is.na(low)
        lst2 <- is.na(upp)
        lst <- lst1 & lst2
        conf.int[lst, ] <- rep(NA, 2)
        ans.all$conf.int <- conf.int
        ans.all$update <- update
        ans.all$profile <- profile.out$profile
        if (quick.ans == 1) {
            cat("\n\nCheck the log-likelihood plot\n")
        }
        ans.all$text.CED <- profile.out$text.CED
        if (update) 
            if (cont) 
                ans.all$show <- f.show.con(ans.all)
            else ans.all$show <- f.show.cat(ans.all)
        if (exists("track")) 
            print("f.CI END")
        return(ans.all)
    })
}



f.profile.all <- function(ans.all, nolog = F, debug = FALSE, display_plots = TRUE) {
    if (exists("track")) 
        print("f.profile.all")
    date.start <- date()
    with(ans.all, {
        if (debug) {
            xlow <- eval(parse(prompt = paste(" give lower bound of CED for plotting profile -------- > ")))
            xupp <- eval(parse(prompt = paste(" give upper bound of CED for plotting profile -------- > ")))
            x.lim <- log10(c(xlow, xupp))
            interrupt <- menu(c("no", "yes"), title = "Do you want to interrupt at each run? \n  ") - 
                1
        }
        crit <- 0.5 * qchisq(conf.lev, 1)
        disturb <- 1.01
        loglik.max <- loglik
        lb.orig <- lb
        ub.orig <- ub
        profile.out <- list()
        tb <- "\t"
        if (dtype == 3) {
        }
        from <- 0
        if (cont) {
            sig2 <- mean(MLE[1:nr.var])
            if (group[1] == 0) {
                if (model.ans == 16) {
                  from <- nr.var + nr.aa + nr.bb + 2
                  until <- from
                  group <- from:until
                }
                else {
                  from <- nr.var + nr.aa + 1
                  if (length(xans) == 1) {
                    until <- nr.var + nr.aa + nr.bb
                    group <- from:until
                  }
                  else {
                    from.rpf <- length(MLE) - nr.dosecol + 2
                    until.rpf <- length(MLE)
                    group <- c(from, from.rpf:until.rpf)
                  }
                }
            }
            if (dtype %in% c(5, 15) && nruns > 0) 
                return(group)
        }
        if (!cont) {
            if (group[1] == 0) {
                from <- nr.aa + 1
                until <- nr.aa + nr.bb
                if (model.type == 2) {
                  if (nr.aa == 1) {
                    until <- until - nr.bb + 1
                  }
                  if (nr.aa > 1 && nr.bb == 1 && ces.ans != 4) {
                    from <- 2
                    until <- nr.aa + 1
                  }
                  if (nr.aa == 1 & nr.bb > 1) {
                    from <- 2
                    until <- nr.bb + 1
                  }
                  if (nr.aa > 1 & nr.bb > 1) {
                    from <- nr.aa + 1
                    until <- nr.aa + nr.bb
                  }
                }
                if (model.ans == 30 & model.type == 1) {
                  from <- nr.aa + nr.bb + 1
                  until <- nr.aa + nr.bb + 1
                }
                if (dtype == 6) {
                  from <- from + 1
                  until <- until + 1
                }
                group <- from:until
            }
        }
        conf.int <- matrix(NA, nrow = length(group), ncol = 2)
        par.nr <- 0
        for (jj in group) {
            par.nr <- par.nr + 1
            CI.NA.low <- FALSE
            CI.NA.upp <- FALSE
            stop <- FALSE
            if (trace) {
                if (model.type == 2 && nr.aa > 1 && nr.bb == 
                  1) 
                  cat("\n\ncalculating C.I.for group", fct1.txt[par.nr], 
                    text.par[jj], "......\n")
                else if (length(xans) > 1) 
                  cat("\ncalculating C.I.for group", fct2.txt[par.nr], 
                    text.par[jj], "......\n")
                else if (nr.bb == 1 || quick.ans == 2) 
                  cat("\nCalculating C.I. ......\n")
                else if (group[1] == 0) 
                  cat("\nCalculating C.I.for group", fct2.txt[par.nr], 
                    text.par[jj], "......\n")
                else cat("\nCalculating C.I.for parameter", jj, 
                  text.par[jj], "......\n")
            }
            CED.low.inf <- FALSE
            CED.upp.inf <- FALSE
            small.step.low <- F
            small.step.upp <- F
            large.step.low <- F
            large.step.upp <- F
            if (trace) {
                cat("\n=========== lower limit =================\n")
                cat(text.par, "loglik", "loglik boundary", sep = tb)
                cat("\n")
                cat(signif(MLE, 3), round(loglik.max, 2), round(loglik.max - 
                  crit, 2), sep = tb)
            }
            lb <- lb.orig
            ub <- ub.orig
            lb[jj] <- MLE[jj]
            if (lb[jj] < 0) 
                nolog <- TRUE
            else nolog <- FALSE
            loglik.low <- rep(loglik.max, 2)
            CED.low <- rep(MLE[jj], 2)
            start <- disturb * MLE
            if (!debug) 
                if (trace.plt && !(dtype %in% 2:3)) {
                  ans.all.plt <- ans.all
                  ans.all.plt$color <- rep(1, length(x))
                  if (from != 0) 
                    ans.all.plt$color[par.nr] <- 3
                  if (display_plots) {
                    f.graph.window(2, nr.gr = nr.gr)
                  }
                  if (!cont & dtype != 3) {
                    f.plot.frq(ans.all.plt)
                  }
                  else {
                    if (ans.all$plt.mns == 1) {
                      ans.all.plt$plt.mns <- 3
                      ans.all.plt$cex.2 <- 2.5
                    }
                    f.plot.con(ans.all.plt)
                  }
                  title(paste("\n\n", text.par[jj]))
                }
            if (debug) {
                f.graph.window(1)
                print("haha")
                y.lim <- c((loglik.max - 10 * crit), loglik.max)
                if (nolog) 
                  xxx <- CED.low
                else xxx <- log10(CED.low)
                print("haha")
                print(x.lim)
                print(xxx)
                print(loglik.low)
                plot(xxx, loglik.low, xlim = x.lim, ylim = y.lim, 
                  ylab = "log(likelihood)", col = 1, pch = 16)
                lines(x.lim, rep(loglik.max - crit, 2), col = 1, 
                  lty = 2)
                h.pr()
            }
            if (cont && dtype %in% c(1, 5, 10)) 
                step.start <- 1.06 + 0.11 * sig2
            else step.start <- 1.08
            step <- step.start
            run <- 0
            max.runs <- 1000
            while (stop == F) {
                run <- run + 1
                if (!nolog) 
                  lb[jj] <- lb[jj]/step
                if (nolog) 
                  lb[jj] <- lb[jj] * step
                ub[jj] <- lb[jj]
                start[jj] <- lb[jj] * disturb
                if (0) {
                  step.tmp <- step
                  while (lb[jj] <= lb.orig[jj]) {
                    if (trace) {
                      cat("\nf.profile.all: value of parameter exceeds lower constraint (", 
                        lb.orig[jj], ")\n")
                      cat("attempting to avoid this by using smaller step\n")
                      cat("it may be needed to refit the model with other constraints for the target parameter \n")
                    }
                    step.tmp <- sqrt(step.tmp)
                    lb[jj] <- lb.orig[jj] * step.tmp
                    ub[jj] <- lb[jj]
                    start[jj] <- lb[jj]
                    step <- step.tmp
                  }
                }
                start[start < lb] <- lb[start < lb]
                start[start > ub] <- ub[start > ub]
                ans.all$par.start <- start
                ans.all$lb <- lb
                ans.all$ub <- ub
                if (cont) 
                  loglik.try <- -f.lik.con(start, x, y, dtype, 
                    fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                    nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                    ttt = ttt, trace = F, fct4 = fct4, fct5 = fct5, 
                    cens.up = cens.up, par.tmp = NA, increase = increase, 
                    x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                    ans.m6.sd = ans.m6.sd, Mx = Mx, x1 = x1, 
                    x2 = x2, cc.inf = cc.inf, incr.decr.no = incr.decr.no)
                if (!cont) 
                  loglik.try <- -f.lik.cat(start, x, y, kk, nn, 
                    dtype, fct1, fct2, nrp, nth, nr.aa, nr.bb, 
                    model.ans, model.type, th.par = th.par, ttt = ttt, 
                    twice = twice, cens.up = cens.up, fct3 = fct3, 
                    ces.ans = ces.ans, CES = CES, CES.cat = CES.cat, 
                    decr.zz = decr.zz, kk.tot = kk.tot, cc.inf = cc.inf, 
                    nn.tot = nn.tot, ref.lev = ref.lev, Mx = Mx, 
                    x1 = x1, x2 = x2, output = TRUE)
                if (!is.finite(loglik.try) || loglik.try < loglik.max - 
                  1e+08) {
                  lb[jj] <- lb[jj] * step
                  step <- sqrt(step)
                  if (trace) 
                    cat("\ntemporary fitting problem, trying to solve this\n")
                  start <- start * runif(length(start), 0.9, 
                    1.1)
                }
                else {
                  nlminb.out <- f.nlminb(ans.all, tmp.quick = T)
                  MLE.current <- nlminb.out$MLE
                  loglik.current <- nlminb.out$loglik
                  start <- MLE.current * disturb
                  if (trace) {
                    cat("\n")
                    cat(signif(MLE.current, 3), round(loglik.current, 
                      2), round(loglik.max - crit, 2), sep = tb)
                  }
                  if (debug) {
                    if (nolog) 
                      xxx <- MLE.current[jj]
                    else xxx <- log10(MLE.current[jj])
                    points(xxx, loglik.current)
                    cat("\n")
                    print("MLE current")
                    print(MLE.current[jj])
                    print("loglik current")
                    print(loglik.current)
                    if (interrupt) 
                      f.press.key.to.continue()
                  }
                  if (!debug && trace.plt && (loglik.current > 
                    (loglik.max - crit))) {
                    if (!cont & !(dtype %in% 2:3)) {
                      par.lst <- f.split.par(MLE.current, nrp, 
                        nth, dtype, fct3)
                      ans.all.plt$regr.par <- par.lst$regr.par
                      f.lines.frq(ans.all.plt)
                    }
                    else if (cont) {
                      if (ans.m6.sd != 2) 
                        ans.all.plt$regr.par <- MLE.current[-c(1:nr.var)]
                      f.lines.con(ans.all.plt)
                      if (0) {
                        ans.all.plt$regr.par.matr <- f.pars(ans.all.plt)$regr.par.matr
                        ans.all.plt$CED <- f.ced.con(ans.all.plt)
                        f.cedlines.con(ans.all.plt)
                      }
                    }
                  }
                  if (loglik.current > loglik.max + 0.03) {
                    ans.all$par.start <- MLE.current * disturb
                    if (trace) 
                      cat("\nlocal optimum found....\n")
                    profile.out$MLE.new <- nlminb.out$MLE
                    profile.out$loglik <- loglik.current
                    profile.out$conf.int <- conf.int
                    profile.out$converged <- nlminb.out$converged
                    return(profile.out)
                  }
                  loglik.low <- c(loglik.current, loglik.low)
                  CED.low <- c(lb[jj], CED.low)
                  if ((CED.low[1] < CED.low[2]) && loglik.low[1] > 
                    loglik.low[2]) {
                    CED.low <- CED.low[-2]
                    loglik.low <- loglik.low[-2]
                  }
                  if (((loglik.max - loglik.current) > crit) && 
                    (length(unique(loglik.low)) > 6)) {
                    stop <- TRUE
                  }
                  if (step < 1.00001) {
                    if (trace) 
                      cat("\nstep size too small, calculations stopped\n")
                    stop <- T
                    small.step.low <- T
                  }
                  if (step > 1e+10) {
                    if (trace) 
                      cat("\nstep size too large, calculations stopped\n")
                    stop <- T
                    large.step.low <- T
                  }
                  if (abs(CED.low[1]) < 1e-20) {
                    if (trace) 
                      cat("\ncurrent parameter value is too small:", 
                        CED.low[1], ", evaluating likelihood profile is stopped\n")
                    CED.low.inf <- TRUE
                    stop <- T
                  }
                  if (run > max.runs) {
                    if (trace) 
                      cat("\nATTENTION: Maximum number of runs reached in gauching profile\n")
                    stop <- T
                  }
                  if (!stop) {
                    if (loglik.low[2] - loglik.low[1] > 0.25 * 
                      crit) {
                      if (!nolog) 
                        lb[jj] <- lb[jj] * step
                      if (nolog) 
                        lb[jj] <- lb[jj]/step
                      step <- step^0.4
                      CED.low[1] <- NA
                      if (trace) 
                        cat("\nLarge change in loglik, step decreased:", 
                          signif(step, 4), "\n")
                    }
                    else if (loglik.low[2] - loglik.low[1] < 
                      0.1 * crit) {
                      if (loglik.low[2] >= loglik.low[1]) {
                        step <- step^1.26
                        if (!is.na(loglik.low[3])) 
                          if (loglik.low[3] >= loglik.low[1]) 
                            step <- step^1.26
                        if (!is.na(loglik.low[4])) 
                          if (loglik.low[4] >= loglik.low[1]) 
                            step <- step^1.26
                        if (!is.na(loglik.low[5])) 
                          if (loglik.low[5] >= loglik.low[1]) 
                            step <- step^1.26
                        if (trace) 
                          cat("\nSmall change in loglik, step increased:", 
                            signif(step, 4), "\n")
                      }
                    }
                    else if (loglik.low[1] < loglik.max - crit) {
                      if (trace) 
                        cat("\ncritical log-likelihood reached too early, jump back")
                      if (!nolog) 
                        lb[jj] <- MLE.current[jj]/sqrt(step)
                      if (nolog) 
                        lb[jj] <- MLE.current[jj] * sqrt(step)
                      CED.low[1] <- NA
                    }
                  }
                }
                loglik.low <- loglik.low[order(CED.low)]
                CED.low <- sort(CED.low)
            }
            if (trace) {
                cat("\n\n=========== upper limit =================\n")
                cat(text.par, "loglik", "loglik boundary", sep = tb)
                cat("\n")
                cat(signif(MLE, 3), round(loglik.max, 2), round(loglik.max - 
                  crit, 2), sep = tb)
            }
            lb <- lb.orig
            ub <- ub.orig
            lb[jj] <- MLE[jj]
            loglik.upp <- rep(loglik.max, 2)
            CED.upp <- rep(MLE[jj], 2)
            start <- disturb * MLE
            if (cont && dtype %in% c(1, 5, 10)) 
                step.start <- 1.06 + 0.11 * sig2
            else step.start <- 1.08
            step <- step.start
            run <- 0
            stop <- F
            while (stop == F) {
                run <- run + 1
                if (!nolog) 
                  lb[jj] <- lb[jj] * step
                if (nolog) 
                  lb[jj] <- lb[jj]/step
                ub[jj] <- lb[jj]
                start[jj] <- lb[jj]
                if (0) {
                  step.tmp <- step
                  while (lb[jj] >= ub.orig[jj]) {
                    if (trace) {
                      cat("\nf.profile.all: value of parameter exceeds upper constraint (", 
                        lb.orig[jj], ")\n")
                      cat("attempting to avoid this by using smaller step\n")
                      cat("it may be needed to refit the model with other constraints for the target parameter \n")
                    }
                    step.tmp <- sqrt(step.tmp)
                    lb[jj] <- lb[jj]/step.tmp
                    ub[jj] <- lb[jj]
                    step <- step.tmp
                  }
                }
                start[start < lb] <- lb[start < lb]
                start[start > ub] <- ub[start > ub]
                ans.all$par.start <- start
                ans.all$lb <- lb
                ans.all$ub <- ub
                if (cont) 
                  loglik.try <- -f.lik.con(start, x, y, dtype, 
                    fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                    nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                    ttt = ttt, trace = F, fct4 = fct4, fct5 = fct5, 
                    cens.up = cens.up, par.tmp = NA, increase = increase, 
                    x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                    ans.m6.sd = ans.m6.sd, Mx = Mx, x1 = x1, 
                    x2 = x2, cc.inf = cc.inf, incr.decr.no = incr.decr.no)
                if (!cont) 
                  loglik.try <- -f.lik.cat(start, x, y, kk, nn, 
                    dtype, fct1, fct2, nrp, nth, nr.aa, nr.bb, 
                    model.ans, model.type, th.par = th.par, ttt = ttt, 
                    twice = twice, cens.up = cens.up, fct3 = fct3, 
                    ces.ans = ces.ans, CES = CES, CES.cat = CES.cat, 
                    decr.zz = decr.zz, kk.tot = kk.tot, cc.inf = cc.inf, 
                    nn.tot = nn.tot, ref.lev = ref.lev, Mx = Mx, 
                    x1 = x1, x2 = x2, output = TRUE)
                if (!is.finite(loglik.try) || loglik.try < loglik.max - 
                  1e+08) {
                  if (!nolog) 
                    lb[jj] <- lb[jj]/step
                  if (nolog) 
                    lb[jj] <- lb[jj] * step
                  step <- sqrt(step)
                  if (trace) 
                    cat("\ntemporary fitting problem\n")
                  start <- start * runif(length(start), 0.9, 
                    1.1)
                }
                else {
                  nlminb.out <- f.nlminb(ans.all, tmp.quick = T)
                  MLE.current <- nlminb.out$MLE
                  loglik.current <- nlminb.out$loglik
                  start <- MLE.current * disturb
                  if (trace) {
                    cat("\n")
                    cat(signif(MLE.current, 3), round(loglik.current, 
                      2), round(loglik.max - crit, 2), sep = tb)
                  }
                  if (debug) {
                    if (nolog) 
                      xxx <- MLE.current[jj]
                    else xxx <- log10(MLE.current[jj])
                    points(xxx, loglik.current)
                    cat("\n")
                    print(MLE.current[jj])
                    print(loglik.current)
                    if (interrupt) 
                      f.press.key.to.continue()
                  }
                  if (!debug && trace.plt & (loglik.current > 
                    (loglik.max - crit))) {
                    if (!cont && !(dtype %in% 2:3)) {
                      par.lst <- f.split.par(MLE.current, nrp, 
                        nth, dtype, fct3)
                      ans.all.plt$regr.par <- par.lst$regr.par
                      f.lines.frq(ans.all.plt)
                    }
                    else if (cont) {
                      if (ans.m6.sd != 2) 
                        ans.all.plt$regr.par <- MLE.current[-c(1:nr.var)]
                      f.lines.con(ans.all.plt)
                    }
                  }
                  if (loglik.current > loglik.max + 0.03) {
                    ans.all$par.start <- MLE.current * disturb
                    if (trace) 
                      cat("\nlocal optimum found....\n")
                    profile.out$MLE.new <- nlminb.out$MLE
                    profile.out$loglik <- loglik.current
                    profile.out$conf.int <- conf.int
                    profile.out$converged <- nlminb.out$converged
                    return(profile.out)
                  }
                  loglik.upp <- c(loglik.current, loglik.upp)
                  CED.upp <- c(lb[jj], CED.upp)
                  if ((CED.upp[1] > CED.upp[2]) && loglik.upp[1] > 
                    loglik.upp[2]) {
                    CED.upp <- CED.upp[-2]
                    loglik.upp <- loglik.upp[-2]
                  }
                  if ((loglik.max - loglik.current > crit) && 
                    (length(unique(loglik.upp)) > 6)) 
                    stop <- TRUE
                  if (step < 1.00001) {
                    stop <- T
                    small.step.upp <- T
                    if (trace) 
                      cat("\nstep size too small, calculations stopped\n")
                  }
                  if (step > 1e+10) {
                    stop <- T
                    large.step.upp <- T
                    if (trace) 
                      cat("\nstep size too large, calculations stopped\n")
                  }
                  if (abs(CED.upp[1]) > 1e+20) {
                    if (trace) 
                      cat("\ncurrent parameter value is too large:", 
                        CED.upp[1], "\n")
                    CED.upp.inf <- TRUE
                    stop <- T
                  }
                  if (run > max.runs) {
                    if (trace) 
                      cat("\nATTENTION: Maximum number of runs reached in establishing profile\n")
                    stop <- T
                  }
                  if (!stop) {
                    if (loglik.upp[2] - loglik.upp[1] > 0.25 * 
                      crit) {
                      if (!nolog) 
                        lb[jj] <- lb[jj]/step
                      if (nolog) 
                        lb[jj] <- lb[jj] * step
                      step <- step^0.4
                      CED.upp[1] <- NA
                      if (trace) 
                        cat("\nLarge change in loglik, step decreased:", 
                          signif(step, 4), "\n")
                    }
                    else if (loglik.upp[2] - loglik.upp[1] < 
                      0.1 * crit) {
                      if (loglik.upp[2] >= loglik.upp[1]) {
                        step <- step^1.26
                        if (!is.na(loglik.upp[3])) 
                          if (loglik.upp[3] >= loglik.upp[1]) 
                            step <- step^1.26
                        if (!is.na(loglik.upp[4])) 
                          if (loglik.upp[4] >= loglik.upp[1]) 
                            step <- step^1.26
                        if (!is.na(loglik.upp[5])) 
                          if (loglik.upp[5] >= loglik.upp[1]) 
                            step <- step^1.26
                        if (trace) 
                          cat("\nSmall change in loglik, step increased:", 
                            signif(step, 4), "\n")
                      }
                    }
                    else if (loglik.upp[1] < loglik.max - crit) {
                      if (trace) 
                        cat("\ncritical log-likelihood reached too early, jump back")
                      if (!nolog) 
                        lb[jj] <- MLE.current[jj] * sqrt(step)
                      if (nolog) 
                        lb[jj] <- MLE.current[jj]/sqrt(step)
                      CED.upp[1] <- NA
                    }
                    if (0) 
                      if (model.ans %in% c(14, 15, 24, 25)) {
                        cc <- MLE.current[nr.var + nr.aa + nr.bb + 
                          1]
                        if (increase == 1) 
                          if (cc <= 1 + CES) {
                            print("f.profile.all:   c does not reach CES")
                            print(MLE.current)
                            print(cc)
                            stop <- T
                            CED.upp.inf <- TRUE
                          }
                        if (increase == -1) 
                          if (cc >= 1 - abs(CES)) {
                            print("f.profile.all decrease:   c does not reach CES")
                            print(MLE.current)
                            print(cc)
                            print(1 - abs(CES))
                            stop <- T
                            CED.upp.inf <- TRUE
                          }
                      }
                  }
                }
                if (is.na(CED.upp[1])) {
                  CED.upp <- CED.upp[-1]
                  loglik.upp <- loglik.upp[-1]
                }
                loglik.upp <- loglik.upp[rev(order(CED.upp))]
                CED.upp <- rev(sort(CED.upp))
            }
            if (debug) 
                f.press.key.to.continue()
            lst <- loglik.low < (loglik.max - crit)
            if (sum(lst) > 0) {
                first.crit.ll <- max(loglik.low[lst])
                first.crit.ced <- CED.low[loglik.low == first.crit.ll]
                first.crit.ll <- rep(first.crit.ll, length(first.crit.ced))
                loglik.low <- c(loglik.low[!lst], first.crit.ll)
                CED.low <- c(CED.low[!lst], first.crit.ced)
                loglik.low <- loglik.low[order(CED.low)]
                CED.low <- sort(CED.low)
            }
            lst <- loglik.upp < (loglik.max - crit)
            if (length(loglik[lst]) > 0) {
                first.crit.ll <- max(loglik.upp[lst])
                first.crit.ced <- CED.upp[loglik.upp == first.crit.ll]
                first.crit.ll <- rep(first.crit.ll, length(first.crit.ced))
                loglik.upp <- c(loglik.upp[!lst], first.crit.ll)
                CED.upp <- c(CED.upp[!lst], first.crit.ced)
                loglik.upp <- loglik.upp[order(CED.upp)]
                CED.upp <- sort(CED.upp)
            }
            if (!nolog) 
                CED.low <- log10(CED.low)
            CED.low <- CED.low[-length(CED.low)]
            loglik.low <- loglik.low[-length(loglik.low)]
            if (any(diff(loglik.low) < 0)) 
                spline.low <- spline(CED.low, loglik.low)
            else spline.low <- spline(CED.low, loglik.low, method = "hyman")
            if (min(spline.low$x) > min(CED.low)) {
                spline.low$x <- c(min(CED.low), spline.low$x)
                spline.low$y <- c(loglik.low[1], spline.low$y)
            }
            if (length(unique(spline.low$y)) == 1) {
                if (trace) 
                  cat("\n\nonly one point in lower spline\n")
                ci.low <- list(y = -Inf)
                conf.int[par.nr, 1] <- ci.low$y
                CI.NA.low <- TRUE
            }
            else {
                ci.low <- approx(spline.low$y, spline.low$x, 
                  xout = loglik.max - crit)
                conf.int[par.nr, 1] <- ci.low$y
            }
            if (small.step.low) {
                conf.int[par.nr, 1] <- CED.low[1]
            }
            if (large.step.low) 
                conf.int[par.nr, 1] <- -Inf
            if (loglik.low[1] > loglik.max - crit) 
                conf.int[par.nr, 1] <- -Inf
            if (!nolog) 
                CED.upp <- log10(CED.upp)
            if (max(abs(diff(loglik.upp))) > 0.01) {
                CED.upp <- CED.upp[-1]
                loglik.upp <- loglik.upp[-1]
                if (any(diff(loglik.upp) < 0)) 
                  spline.upp <- spline(CED.upp, loglik.upp)
                else spline.upp <- spline(CED.upp, loglik.upp, 
                  method = "hyman")
                if (max(spline.upp$x) < max(CED.upp)) {
                  spline.upp$x <- c(spline.upp$x, max(CED.upp))
                  spline.upp$y <- c(spline.upp$y, loglik.upp[length(loglik.upp)])
                }
                if (length(unique(spline.upp$y)) == 1) {
                  if (trace) 
                    cat("\n\nonly one point in upper spline\n")
                  ci.upp <- list(y = Inf)
                  conf.int[par.nr, 2] <- ci.upp$y
                  CI.NA.upp <- TRUE
                }
                else {
                  ci.upp <- approx(spline.upp$y, spline.upp$x, 
                    xout = loglik.max - crit)
                  conf.int[par.nr, 2] <- ci.upp$y
                }
                if (small.step.upp) 
                  conf.int[par.nr, 2] <- CED.upp[1]
                if (large.step.upp) 
                  conf.int[par.nr, 2] <- Inf
                if (loglik.upp[length(loglik.upp)] > loglik.max - 
                  crit) {
                  if (!small.step.upp) {
                    conf.int[par.nr, 2] <- Inf
                  }
                  else conf.int[par.nr, 2] <- NA
                }
            }
            else {
                spline.upp <- list()
                spline.upp$x <- CED.upp
                spline.upp$y <- rep(loglik.upp[1], length(CED.upp))
                ci.upp <- list(y = max(CED.upp))
            }
            if (CED.low.inf) 
                conf.int[par.nr, 1] <- -Inf
            if (CED.upp.inf) 
                conf.int[par.nr, 2] <- Inf
            CED.vec <- c(CED.low, CED.upp)
            loglik.vec <- c(loglik.low, loglik.upp)
            if (!(CI.NA.low && CI.NA.upp)) 
                if (trace.plt) {
                  if (dtype %in% 2:3)
                    if (display_plots) {
                      f.graph.window(1, nr.gr = nr.gr)
                    } 
                    y.low <- min(spline.low$y, spline.upp$y)
                  y.upp <- max(spline.low$y, spline.upp$y)
                  y.lim <- c(min(y.low, loglik.max - crit), y.upp)
                  if (nolog) 
                    x.lab = text.par[jj]
                  else x.lab = paste("log10(", text.par[jj], 
                    ")")
                  plot(CED.vec, loglik.vec, xlab = x.lab, ylim = y.lim, 
                    ylab = "log(likelihood)", col = 1, pch = 16)
                  title("calculation of confidence interval")
                  if (nolog) 
                    MLE.tmp <- MLE[jj]
                  else MLE.tmp <- log10(MLE[jj])
                  points(MLE.tmp, loglik.max, col = 8, pch = 16, 
                    cex = 2.5)
                  lines(spline.low$x, spline.low$y, col = 6)
                  lines(spline.upp$x, spline.upp$y, col = 6)
                  segments(min(CED.vec), loglik.max - crit, max(CED.vec), 
                    loglik.max - crit, col = 1, lty = 2)
                  lines(rep(ci.low$y, 2), c(min(loglik.vec), 
                    loglik.max), lty = 2)
                  lines(rep(ci.upp$y, 2), c(min(loglik.vec), 
                    loglik.max), lty = 2)
                }
            if (max(abs(diff(loglik.low))) < 0.01) 
                conf.int[par.nr, 1] <- -Inf
            if (max(abs(diff(loglik.upp))) < 0.01) {
                conf.int[par.nr, 2] <- Inf
            }
            if (CI.NA.low) {
                if (trace) 
                  cat("\nlower conf. bound could not be established, log-likelihood profile did not change\n")
            }
            if (CI.NA.upp) {
                if (trace) 
                  cat("\nupper conf. bound could not be established, log-likelihood profile did not change\n")
            }
        }
        if (length(conf.int[, 1]) > 1) {
            text.CED <- text.par[group]
            dimnames(conf.int)[[1]] <- text.CED
        }
        else text.CED <- rep(NA, length(conf.int[, 1]))
        ans.all$text.CED <- text.CED
        if (nolog) 
            conf.int <- signif(conf.int, 3)
        else conf.int <- signif(10^conf.int, 3)
        profile.out <- list(conf.int = conf.int, MLE = MLE, loglik = nlminb.out$loglik, 
            profile = cbind(CED.vec, loglik.vec), text.CED = text.CED)
        if (quick.ans == 1) {
            if (trace) 
                cat("\nstarted at:", date.start)
            if (trace) 
                cat("\n")
            if (trace) 
                cat(date())
        }
        profile.out$CED.low <- CED.low
        profile.out$loglik.low <- loglik.low
        profile.out$CED.upp <- CED.upp
        profile.out$loglik.upp <- loglik.upp
        ans.all$boot <- FALSE
        if (exists("track")) 
            print("end of profile.all")
        return(profile.out)
    })
}



f.plot.gui <- function(ans.all, HTML = FALSE, model.summ = TRUE, display_plots = TRUE, results_env = NULL, output_type = NULL, filename = NULL) {
    if (exists("track")) 
        print("f.plot.gui")
    WAPP <- ans.all$WAPP
    if (ans.all$dtype == 3) 
        ans.all$displ.no <- 0
    if (ans.all$cont || ans.all$dtype == 3) {
        if (ans.all$plot.type %in% c(9, 10)) {
            if (ans.all$dtype %in% c(10, 15, 250, 260)) {
                cat("\n residual plots not possible for summary data\n")
                return(invisible())
            }
            if (ans.all$model.fam == 0) {
                cat("\n\nATTENTION: first select a single model when plotting residuals\n\n")
                return(invisible())
            }
        }
        if (length(ans.all$xans) > 1) {
            ans.all$plt.mns <- 3
            ans.all$CI.plt <- TRUE
        }
        if (ans.all$model.fam == 0 && length(ans.all$SINGMOD) == 
            0) {
            if (!HTML) {
                if (ans.all$nr.models == 1)
                  f.graph.window(1)
                else {
                  if (WAPP & !ans.all$svg.plots) 
                    .gw.size <- c(800, 400)
                  else if (WAPP & ans.all$svg.plots) 
                    .gw.size <- c(8, 4)
                  else .gw.size <- c(7, 3)
                  name.wapp <- paste("a", ans.all$yans, "exphill", 
                    sep = "")
                  if (display_plots) {
                    f.create.graphwin(.gw.size[1], .gw.size[2], 
                    WAPP = WAPP, title = "", name.wapp = name.wapp, 
                    plotprefix = ans.all$plotprefix, svg.plots = ans.all$svg.plots,
                    output_type = output_type, filename = filename)
                  }
                  par(mfcol = c(1, 2))
                  par(mar = c(5.5, 4.8, 4, 12))
                  par(cex.main = 1.2)
                  par(cex.sub = 1)
                  par(cex.lab = 1.5)
                  par(cex = 0.6)
                }
                y.pos <- 0
                assign(".ypos", y.pos, immediate = T, pos = 1)
            }
            ans.all.plt <- ans.all
            ans.all.plt$heading <- ""
            if (ans.all$EXP$trend == FALSE) 
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$NULL.model)
            else ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$EXP)
            if (model.summ) 
                if (ans.all$dtype == 3) 
                  ans.all.plt$show <- f.show.cat(ans.all.plt)
                else ans.all.plt$show <- f.show.con(ans.all.plt)
            ans.all <- f.plot.all(ans.all.plt, new.window = F, display_plots = display_plots)
            if (ans.all$dtype == 3) 
                title(main = "Exponential model", font.main = 2.5)
            if (length(ans.all$xans) > 1 && max(ans.all$fct1) > 
                1) 
                title(main = "\n\nExponential model", font.main = 2.5)
            if (ans.all$nr.models > 1) {
                ans.all.plt <- ans.all
                if (ans.all$HILL$trend == FALSE) 
                  ans.all.plt <- f.move.sublist(ans.all.plt, 
                    ans.all$NULL.model)
                else ans.all.plt <- f.move.sublist(ans.all.plt, 
                  ans.all$HILL)
                y.pos <- 95
                assign(".ypos", y.pos, immediate = T, pos = 1)
                ans.all.plt$yleg <- ""
                if (model.summ) 
                  if (ans.all$dtype == 3) 
                    ans.all.plt$show <- f.show.cat(ans.all.plt)
                  else ans.all.plt$show <- f.show.con(ans.all.plt)
                f.plot.all(ans.all.plt, new.window = F, display_plots = display_plots)
                if (ans.all$dtype == 3) 
                  title(main = "Hill model", font.main = 2)
                if (length(ans.all$xans) > 1 && max(ans.all$fct1) > 
                  1) 
                  title(main = "\n\nHill model", font.main = 2)
            }
            if (ans.all$nr.models > 2) {
                if (WAPP) 
                  dev.off()
                name.wapp <- paste("a", ans.all$yans, "invlogn", 
                  sep = "")
                if (display_plots) {
                  f.create.graphwin(.gw.size[1], .gw.size[2], WAPP = WAPP, 
                    title = "", name.wapp = name.wapp, plotprefix = ans.all$plotprefix, 
                    svg.plots = ans.all$svg.plots,
                    output_type = output_type, filename = filename)
                  par(mfcol = c(1, 2))
                  par(mar = c(5.5, 4.8, 4, 12))
                  par(cex.main = 1.2)
                  par(cex.sub = 1)
                  par(cex.lab = 1.5)
                  par(cex = 0.6)
                }
            }
            if (ans.all$nr.models > 2) {
                ans.all.plt <- ans.all
                if (ans.all$INVEXP$trend == FALSE) 
                  ans.all.plt <- f.move.sublist(ans.all.plt, 
                    ans.all$NULL.model)
                else ans.all.plt <- f.move.sublist(ans.all.plt, 
                  ans.all$INVEXP)
                y.pos <- 95
                assign(".ypos", y.pos, immediate = T, pos = 1)
                ans.all.plt$yleg <- ""
                if (model.summ) 
                  if (ans.all$dtype == 3) 
                    ans.all.plt$show <- f.show.cat(ans.all.plt)
                  else ans.all.plt$show <- f.show.con(ans.all.plt)
                f.plot.all(ans.all.plt, new.window = F, display_plots = display_plots)
                if (ans.all$dtype == 3) 
                  title(main = "inverse expon model", font.main = 2)
                if (length(ans.all$xans) > 1 && max(ans.all$fct1) > 
                  1) 
                  title(main = "\n\ninverse expon model", font.main = 2)
            }
            if (ans.all$nr.models > 3) {
                ans.all.plt <- ans.all
                if (ans.all$LOGN$trend == FALSE) 
                  ans.all.plt <- f.move.sublist(ans.all.plt, 
                    ans.all$NULL.model)
                else ans.all.plt <- f.move.sublist(ans.all.plt, 
                  ans.all$LOGN)
                y.pos <- 95
                assign(".ypos", y.pos, immediate = T, pos = 1)
                ans.all.plt$yleg <- ""
                if (model.summ) 
                  if (ans.all$dtype == 3) 
                    ans.all.plt$show <- f.show.cat(ans.all.plt)
                  else ans.all.plt$show <- f.show.con(ans.all.plt)
                f.plot.all(ans.all.plt, new.window = F, display_plots = display_plots)
                if (ans.all$dtype == 3) 
                  title(main = "LN model", font.main = 2)
                if (length(ans.all$xans) > 1 && max(ans.all$fct1) > 
                  1) 
                  title(main = "\n\nLN model", font.main = 2)
            }
        }
        if (ans.all$model.fam > 0 || length(ans.all$SINGMOD) > 
            0) {
            if (HTML || WAPP) {
                new.window <- F
                name.wapp <- paste("a", ans.all$yans, "singmod")
                f.graph.window(1, WAPP = WAPP, name.wapp = name.wapp, 
                  plotprefix = ans.all$plotprefix, svg.plots = ans.all$svg.plots)
            }
            else {
                new.window <- T
                y.pos <- 0
                assign(".ypos", y.pos, immediate = T, pos = 1)
            }
            ans.all.plt <- ans.all
            if (length(ans.all$SINGMOD) > 0) 
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$SINGMOD)
            else switch(ans.all$model.fam, {
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$EXP)
                title.tmp <- "Exponential model"
            }, {
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$HILL)
                title.tmp <- "Hill model"
            }, {
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$INVEXP)
                title.tmp <- "Inverse exponential model"
            }, {
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$LOGN)
                title.tmp <- "Lognormal DR model"
            })
            if (0) 
                if (ans.all.plt$trend == FALSE) 
                  ans.all.plt <- f.move.sublist(ans.all.plt, 
                    ans.all$NULL.model)
            ans.all.plt$heading <- ""
            if (model.summ) 
                if (ans.all$dtype == 3) 
                  ans.all.plt$show <- f.show.cat(ans.all.plt)
                else ans.all.plt$show <- f.show.con(ans.all.plt)
            f.plot.all(ans.all.plt, new.window = new.window, display_plots = display_plots)
            if (ans.all$dtype == 3) 
                title(main = title.tmp, font.main = 2.5)
        }
    }
    if (!ans.all$cont && ans.all$dtype != 3) {
        if (ans.all$gui && !ans.all$WAPP) {
            if (ans.all$model.fam > 0) {
                ans.all$model.ans <- ans.all$model.fam
            }
        }
        if (ans.all$model.fam == 0) {
            if (!HTML) {
                name.wapp <- paste("a", ans.all$yans, "setofmod", 
                  sep = "")
                f.graph.window(9, WAPP = WAPP, title = "quantal models", 
                  name.wapp = name.wapp, plotprefix = ans.all$plotprefix, 
                  svg.plots = ans.all$svg.plots)
            }
            count <- 0
            nr.models <- length(ans.all$model.list) - 2
            for (ii in 1:nr.models) {
                ans.all.plt <- ans.all
                name <- ans.all$model.list[ii]
                if (0) {
                  if (ans.all[[name]]$trend) 
                    ans.all.plt <- f.move.sublist(ans.all.plt, 
                      ans.all[[name]])
                  else {
                    ans.all.plt <- f.move.sublist(ans.all.plt, 
                      ans.all$NULL.model)
                    ans.all.plt$modelname <- ans.all[[name]]$modelname
                  }
                }
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all[[name]])
                ans.all.plt$heading <- ""
                if (ans.all$dtype == 2) 
                  ans.all.plt$CI.plt <- FALSE
                count <- count + 1
                xleg <- ans.all.plt$xleg
                yleg <- ans.all.plt$y.leg
                f.plot.all(ans.all.plt, new.window = FALSE, display_plots = display_plots)
                title(main = paste(ans.all.plt$modelname, "--"), 
                  cex = 0.5)
            }
            if (ans.all$LVM) {
                ans.all.plt <- ans.all
                ans.all.plt$model.type <- 2
                if (0) 
                  if (ans.all$EXP$trend) 
                    ans.all.plt <- f.move.sublist(ans.all.plt, 
                      ans.all$EXP)
                  else {
                    ans.all.plt <- f.move.sublist(ans.all.plt, 
                      ans.all$NULL.model)
                    ans.all.plt$modelname <- ans.all[[name]]$modelname
                  }
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$EXP)
                ans.all.plt$heading <- ""
                if (ans.all$dtype == 2) 
                  ans.all.plt$CI.plt <- FALSE
                count <- count + 1
                f.plot.all(ans.all.plt, new.window = FALSE, display_plots = display_plots)
                title(main = ans.all$EXP$modelname, cex = 0.5)
                ans.all.plt <- ans.all
                ans.all.plt$model.type <- 2
                if (1) 
                  ans.all.plt <- f.move.sublist(ans.all.plt, 
                    ans.all$HILL)
                else ans.all.plt <- f.move.sublist(ans.all.plt, 
                  ans.all$NULL.model)
                ans.all.plt$heading <- ""
                if (ans.all$dtype == 2) 
                  ans.all.plt$CI.plt <- FALSE
                count <- count + 1
                f.plot.all(ans.all.plt, new.window = FALSE, display_plots = display_plots)
                title(main = ans.all$HILL$modelname, cex = 0.5)
            }
            if (count < 9) {
                plot(1:3, 1:3, col = 0, yaxt = "n", xlab = "", 
                  ylab = "", type = "n")
                text(1.2, 2.9, ans.all$odt.name, adj = 0)
                text(1.2, 2.5, "x-axis:", adj = 0)
                text(1.2, 2.3, xleg, adj = 0)
                text(1.2, 2.1, "y-axis:", adj = 0)
                text(1.2, 1.9, yleg, adj = 0)
                text(1.2, 1.2, paste("PROAST version", ans.all$PRversion), 
                  adj = 0)
            }
            else mtext(ans.all$xleg, side = 1)
        }
        if (ans.all$model.fam > 0) {
            if (!HTML) 
                name.wapp <- paste("a", ans.all$yans, "singmod")
            f.graph.window(1, WAPP = WAPP, name.wapp = name.wapp, 
                plotprefix = ans.all$plotprefix, svg.plots = ans.all$svg.plots)
            if (ans.all$model.ans %in% 1:9) {
                ii <- ans.all$model.ans
                ans.all.plt <- ans.all
                name <- ans.all$model.list[ii]
                ans.all.plt <- f.move.sublist(ans.all.plt, ans.all[[name]])
                if (ans.all$dtype == 2) 
                  ans.all.plt$CI.plt <- FALSE
                if (model.summ) 
                  ans.all.plt$show <- f.show.cat(ans.all.plt)
                f.plot.all(ans.all.plt, new.window = FALSE, display_plots = display_plots)
            }
        }
    }
    if (ans.all$dtype != 3) 
        if (length(ans.all$gr.txt) > 1) 
            f.explain.marks(ans.all)
    if (WAPP) 
        dev.off()
    if (exists("track")) 
        print("f.plot.gui:  END")
    if (!is.null(results_env)) {
      if (length(ans.all$Vyans) > 1) {
        ans.all$res.name <- ans.all$varnames[ans.all$yans]
        assign(paste(ans.all$res.name, "plot_result"), ans.all, envir = results_env) 
      } else {
        assign("plot_result", ans.all, envir = results_env)
      }
    }
    return(ans.all)
}



f.CI.sel <- function(ans.all, interactive_mode = NULL, results_env = NULL) {
    if (exists("track")) 
        print("f.CI.sel")
    ans.all <- with(ans.all, {
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        nr.var <- max(fct3)
        if (cont) {
            cc.OK <- f.check.cc(ans.all)
            if (!cc.OK) {
                note.tmp <- paste(y.leg, ": parameter cc too close to CES; no BMD CI calculated")
                ans.all$notes <- paste(ans.all$notes, "\n\n", 
                  note.tmp, "\n")
                ans.all$conf.int <- matrix(NA, ncol = 2, nrow = max(fct2))
                return(ans.all)
            }
            from <- nr.var + nr.aa + 1
            if (length(xans) == 1) {
                until <- nr.var + nr.aa + nr.bb
                ans.all$group <- from:until
            }
            else {
                from.rpf <- nr.var + nr.aa + length(MLE) - nr.dosecol + 
                  2
                until.rpf <- nr.var + nr.aa + length(MLE)
                ans.all$group <- c(from, from.rpf:until.rpf)
            }
        }
        if (!cont) {
            CED.lst <- f.ced.cat(ans.all)
            ans.all$CED.matr <- CED.lst$CED.matr
            ans.all$response.matr <- CED.lst$response.matr
            CED <- CED.lst$CED
            ans.all$CED <- CED
            rank.low <- order(CED)[1]
            if (max(fct1) > 1 && max(fct2) == 1) 
                group <- rank.low + 1
            else group <- rank.low + max(fct1)
            if (dtype == 6) 
                group <- group + 1
            if (length(xans) == 1) 
                ans.all$group <- group
            ans.all$rank.low <- rank.low
            if (model.ans > 1) 
                ans.all$sens.lev <- covar.txt[rank.low]
            else ans.all$sens.lev <- NA
        }
        if (dtype %in% c(2, 4, 6)) 
            no.CI <- T
        if (!no.CI) {
            if (dtype %in% c(5, 15)) {
                ans.all$plot.ans <- 1
                ans.all$nrp <- length(regr.par)
                ans.all$CED <- CED
                ans.all <- f.mm7.con(ans.all)
                conf.int <- ans.all$conf.int
            }
            else {
                if (ans.all$loglik == 1e-12) 
                  cat("\n\nno finite log-likelihood found, CI will not be calculated\n\n")
                else {
                  ans.all$trace <- F
                  ans.all$trace.plt <- F
                  if (length(xans) > 1) {
                    group <- nr.var + nr.aa + 1
                    RPF.rank <- (npar - nr.dosecol + 2):npar
                    if (length(xans) > 1) 
                      group <- c(group, RPF.rank)
                    ans.all$group <- group
                  }
                  if (ans.all$output) 
                    cat("\n\n calculating confidence intervals ....\n\n")
                  ans.all <- f.CI(ans.all)
                  conf.int <- ans.all$conf.int
                  if (ans.all$update) {
                    if (cont) {
                      CED <- ans.all$regr.par[(nr.aa + 1):(nr.aa + 
                        nr.bb)]
                    }
                    if (!cont) {
                      CED.lst <- f.ced.cat(ans.all)
                      ans.all$CED.matr <- CED.lst$CED.matr
                      ans.all$response.matr <- CED.lst$response.matr
                      CED <- CED.lst$CED
                      if (nr.bb == 1) 
                        CED <- rep(CED, nr.aa)
                    }
                    ans.all$Vloglik[ans.all$row.sel] <- ans.all$loglik
                    ans.all$Vconverged[ans.all$row.sel] <- ans.all$converged
                    ans.all$regr.par.matr <- f.pars(ans.all)$regr.par.matr
                  }
                }
            }
            if (model.ans > 1) {
                CED <- signif(CED, 4)
                ans.all$conf.int <- ans.all$conf.int
                conf.int <- signif(conf.int, 3)
                CED.uniq <- unique(CED)
                if (cont && ans.all$output) 
                  if (length(CED.uniq) > 1) {
                    if (model.ans == 46) {
                      cat("\n\nthe CED (in orig. units) and the", 
                        100 * conf.lev, "% confidence interval for group", 
                        covar.txt[1], "is: \n", sf.x * CED[1], 
                        "\n", conf.int[1, 1], "\n", conf.int[1, 
                          2], "\n\n")
                      for (qq in 2:max(fct2)) {
                        cat("\n\nthe RPF (in orig. units) and the", 
                          100 * conf.lev, "% confidence interval for group", 
                          covar.txt[qq], "is: \n", sf.x * CED[qq], 
                          "\n", conf.int[qq, 1], "\n", conf.int[qq, 
                            2], "\n\n")
                      }
                    }
                    else {
                      for (qq in 1:max(fct2)) {
                        cat("\n\nthe CED (in orig. units) and the", 
                          100 * conf.lev, "% confidence interval for group", 
                          covar.txt[qq], "is: \n", sf.x * CED[qq], 
                          "\n", conf.int[qq, 1], "\n", conf.int[qq, 
                            2], "\n\n")
                      }
                    }
                  }
                  else cat("\n\nthe CED (in orig. units) and the", 
                    100 * conf.lev, "% confidence interval is: \n", 
                    sf.x * CED[1], "\n", conf.int[1, 1], "\n", 
                    conf.int[1, 2], "\n\n")
            }
            if (length(CED) > 1) 
                ced.table <- data.frame(subgroup = covar.txt, 
                  BMDL = conf.int[, 1], BMDU = conf.int[, 2], 
                  BMD = CED)
            else ced.table <- data.frame(BMDL = conf.int[, 1], 
                BMDU = conf.int[, 2], BMD = CED)
            ans.all$ced.table <- ced.table
        }
        if (0) 
            if (!gui && length(ans.all$MLE) > 10) {
                f.store.results(ans.all, "refit.tmp")
                cat("\n\nNOTE: confidence intervals were stored in refit.tmp\n\n")
            }
        ans.all$CED <- signif(CED, 4)
        if (exists("track")) 
            print("f.CI.sel:  END")
        if (interactive_mode == FALSE) {
          if (length(ans.all$Vyans) > 1) {
            message("Vyans is greater than 1")
            ans.all$res.name <- ans.all$varname[ans.all$yans]
            assign(paste(ans.all$res.name, ans.all$modelname), ans.all, envir = results_env)
          } else {
            ans.all$res.name <- ans.all$varname[ans.all$yans]
            assign(ans.all$modelname, ans.all, envir = results_env)
        }
        }
        return(ans.all)
    })
}



f.CED.all <- function(CED.all, y.leg, exp.out = NA, hill.out = NA, invexp.out = NA, 
    logn.out = NA, TREND) {
    if (exists("track")) 
        print("f.CED.all")
    CED.matr <- CED.all$CED.matr
    endpoints <- CED.all$endpoints
    Mces <- CED.all$Mces
    Msd <- CED.all$Msd
    covar.txt <- CED.all$covar.txt
    Vtrend <- CED.all$Vtrend
    Vtrend <- c(Vtrend, TREND)
    if (length(exp.out) > 0) {
        nr.aa <- max(exp.out$fct1)
        nr.bb <- max(exp.out$fct2)
        CED <- exp.out$CED
        conf.int <- exp.out$conf.int
        if (length(covar.txt) > 2) {
            cat("\nATTENTION: plot of all CED CIs not implemented for more than two subgroups\n")
            return(CED.all)
        }
        else if (length(conf.int[, 1]) == 1 && length(covar.txt) > 
            1) {
            confint.tmp <- conf.int[1, ]
            conf.int <- matrix(ncol = 2, nrow = length(covar.txt))
            for (ii in 1:length(covar.txt)) conf.int[ii, ] <- confint.tmp
        }
        if (!is.na(CED[1])) {
            col <- 1
            CED.tmp <- rep(NA, length(CED.matr[1, ]))
            for (kk in 1:length(covar.txt)) {
                CED.tmp[(col):(col + 1)] <- signif(conf.int[kk, 
                  ], 3)
                col <- col + 2
            }
            CED.matr <- rbind(CED.matr, CED.tmp)
            endpoints <- c(endpoints, y.leg)
        }
        ces.tmp <- c(y.leg, round(exp.out$CES, 3))
        if (!TREND) 
            ces.tmp <- NA
        sd.tmp <- c(y.leg, round(exp.out$pooled.sd, 4))
    }
    if (length(hill.out) > 1) {
        nr.aa <- max(hill.out$fct1)
        nr.bb <- max(hill.out$fct2)
        CED <- hill.out$CED
        conf.int <- hill.out$conf.int
        if (length(conf.int[, 1]) == 1 && length(covar.txt) > 
            1) {
            confint.tmp <- conf.int[1, ]
            conf.int <- matrix(ncol = 2, nrow = length(covar.txt))
            for (ii in 1:length(covar.txt)) conf.int[ii, ] <- confint.tmp
        }
        if (!is.na(CED[1])) {
            col <- 1
            CED.tmp <- rep(NA, length(CED.matr[1, ]))
            for (kk in 1:length(covar.txt)) {
                CED.tmp[(col):(col + 1)] <- signif(conf.int[kk, 
                  ], 3)
                col <- col + 2
            }
            CED.matr <- rbind(CED.matr, CED.tmp)
            endpoints <- c(endpoints, y.leg)
        }
    }
    if (length(invexp.out) > 1) {
        nr.aa <- max(invexp.out$fct1)
        nr.bb <- max(invexp.out$fct2)
        CED <- invexp.out$CED
        conf.int <- invexp.out$conf.int
        if (length(conf.int[, 1]) == 1 && length(covar.txt) > 
            1) {
            confint.tmp <- conf.int[1, ]
            conf.int <- matrix(ncol = 2, nrow = length(covar.txt))
            for (ii in 1:length(covar.txt)) conf.int[ii, ] <- confint.tmp
        }
        if (!is.na(CED[1])) {
            col <- 1
            CED.tmp <- rep(NA, length(CED.matr[1, ]))
            for (kk in 1:length(covar.txt)) {
                CED.tmp[(col):(col + 1)] <- signif(conf.int[kk, 
                  ], 3)
                col <- col + 2
            }
            CED.matr <- rbind(CED.matr, CED.tmp)
            endpoints <- c(endpoints, y.leg)
        }
    }
    if (length(logn.out) > 1) {
        nr.aa <- max(logn.out$fct1)
        nr.bb <- max(logn.out$fct2)
        CED <- logn.out$CED
        conf.int <- logn.out$conf.int
        if (length(conf.int[, 1]) == 1 && length(covar.txt) > 
            1) {
            confint.tmp <- conf.int[1, ]
            conf.int <- matrix(ncol = 2, nrow = length(covar.txt))
            for (ii in 1:length(covar.txt)) conf.int[ii, ] <- confint.tmp
        }
        if (!is.na(CED[1])) {
            col <- 1
            CED.tmp <- rep(NA, length(CED.matr[1, ]))
            for (kk in 1:length(covar.txt)) {
                CED.tmp[(col):(col + 1)] <- signif(conf.int[kk, 
                  ], 3)
                col <- col + 2
            }
            CED.matr <- rbind(CED.matr, CED.tmp)
            endpoints <- c(endpoints, y.leg)
        }
    }
    Mces <- rbind(Mces, ces.tmp)
    Msd <- rbind(Msd, sd.tmp)
    CED.all$endpoints <- endpoints
    CED.all$CED.matr <- CED.matr
    CED.all$Mces <- Mces
    CED.all$Msd <- Msd
    Vtrend <- CED.all$Vtrend
    CED.all$Vtrend <- c(Vtrend, TREND)
    if (exists("track")) 
        print("f.CED.all:   END")
    return(CED.all)
}



f.plot.CED <- function(ans.all,
                       logscale = T,
                       xx.lim = NA,
                       WAPP = FALSE,
                       plotprefix = "", 
                       svg.plots = FALSE,
                       display_plots = TRUE,
                       interactive_mode = TRUE) {
    if (exists("track"))
        print("f.plot.CED")
    data <- ans.all$CED.all
    covar.txt <- data$covar.txt
    nr.lev <- length(covar.txt)
    if (nr.lev == 0)
        nr.lev <- 1
    endp.name <- data$endpoints
    if (ans.all$do.MA)
        endp.name <- unique(endp.name)
    if (0) {
        if (ans.all$NES.ans == 2)
            endp.name.legend <- paste(endp.name, "", data$Mces[, 
                2])
        else endp.name.legend <- endp.name
    }
    endp.name.legend <- paste(unique(endp.name), data$Mces[, 
        2])
    CI <- data$CED.matr
    if (ans.all$do.MA) {
        CI <- ans.all$CED.matr.ma.plt
        CI <- CI[, -1]
    }
    if (is.na(xx.lim[1])) {
        max.x <- max(CI[CI != Inf], na.rm = T)
        zero.dum <- min(CI[CI > 0], na.rm = T)/1.1
        xx.lim <- c(log10(zero.dum), log10(max.x))
    }
    CES <- data$CES
    NES.ans <- data$NES.ans
    if (logscale)
        CI <- log10(CI)
    nr.row <- length(CI[, 1])
    nr.endp <- length(unique(data$endpoints))
    mainjump <- 1
    subjump <- 0.35
    if (ans.all$do.MA)
        max.y <- nr.endp
    else max.y <- nr.endp * (data$nr.models * subjump + mainjump)
    if (max.y > 40)
        cex.mtext <- 0.6
    else cex.mtext <- 0.8
    yy.lim <- c(1, max.y)
    xx.lab <- paste("CED", CES, sep = "-")
    if (NES.ans == 2) 
        xx.lab <- paste("CED (for endpoint-specific CES)")
    if (logscale) 
        xx.lab <- paste("log10-", xx.lab)
    jj <- 1
    for (ii in 1:nr.lev) {
        ci.matr <- CI[, jj:(jj + 1)]
        ci.tmp <- ci.matr[1, 1:2]
        name.wapp <- paste("b", ii, "cedCI", sep = "")
        endp.nr <- 1
        if (display_plots) {
          f.graph.window(1, WAPP = WAPP, title = "", name.wapp = name.wapp, 
              plotprefix = plotprefix, svg.plots = svg.plots)
          plot(1:max.y, 1:max.y, xlim = xx.lim, ylim = yy.lim,
               ylab = "", yaxt = "n", xlab = xx.lab, type = "n")
          if (length(covar.txt) > 1)
              title(main = paste("group", covar.txt[ii]), cex.main = 1.2)
          if (!ans.all$do.MA)
              title(main = paste("\n\n\n", ans.all$nr.models, "models per endpoint"), 
                  cex.main = 0.8)
          if (ans.all$do.MA) 
              title(main = paste("\n\n\n", "CIs based on model averaging"), 
                  cex.main = 0.8)
        }
        yyy.tmp <- max.y
        if (!is.na(ci.tmp[1])) {
            lty.tmp <- 1
            if (ci.tmp[1] == -Inf) {
                lty.tmp <- 2
                ci.tmp[1] <- -100
            }
            if (ci.tmp[2] == Inf) {
                lty.tmp <- 2
                ci.tmp[2] <- 100
            }
            if (display_plots) {
              lines(ci.tmp, rep(yyy.tmp, 2), lty = lty.tmp)
              mtext(endp.name.legend[1], 4, 1, adj = 0, las = 1,
                    at = yyy.tmp, cex = cex.mtext)
            }
        } else {
          if (display_plots) {
            points(xx.lim[2], yyy.tmp, pch = "\"")
            mtext(endp.name.legend[1], 4, 1, adj = 0, las = 1, 
                at = yyy.tmp, cex = cex.mtext, col = 2)
          }
        }
        first.kk <- 2
        for (kk in first.kk:nr.row) {
            ci.tmp <- ci.matr[kk, 1:2]
            if (endp.name[kk] == endp.name[kk - 1])
                yyy.tmp <- yyy.tmp - subjump
            else {
                endp.nr <- endp.nr + 1
                yyy.tmp <- yyy.tmp - 1
                if (!is.na(ci.tmp[1]))
                  color.tmp <- 1
                else {
                  color.tmp <- 2
                }
                if (display_plots) {
                  mtext(endp.name.legend[endp.nr], 4, 1, adj = 0,
                        las = 1, at = yyy.tmp, cex = cex.mtext,
                        col = color.tmp)
                }
            }
            if (!is.na(ci.tmp[1])) {
                lty.tmp <- 1
                if (ci.tmp[1] == -Inf) {
                  lty.tmp <- 2
                  ci.tmp[1] <- -100
                }
                if (ci.tmp[2] == Inf) {
                  lty.tmp <- 2
                  ci.tmp[2] <- 100
                }
                if (display_plots) {
                  lines(ci.tmp, rep(yyy.tmp, 2), lty = lty.tmp)
                }
            }
            else {
                if (display_plots) {
                  points(xx.lim[2], yyy.tmp, pch = "\"")
                }
            }
        }
        jj <- jj + 2
        if (!WAPP) {
            cat("\nCED-CI plot created for subgroup", ii, "\n")
            if (ii < nr.lev)
                if (interactive_mode == TRUE) {
                  f.press.key.to.continue()
                }
        }
        if (WAPP)
            dev.off()
    }
    if (exists(" track"))
        print("f.plot.CED: END")
}



parse_PROAST_output <- function(result) {
  result <- result[!grepl("plot_result", names(result))]
  selected_models <- c()
  CES <- c()
  CED <- c()
  CEDL <- c()
  CEDU <- c()
  AIC <- c()
  log_likelihood <- c()
  var <- c()
  a <- c()
  d <- c()
  covariates <- c()
  ma_rows <- list()
  model_weights <- list()
  response <- c()

  # Loop through the list excluding 'model_averaging'
  for (i in seq_along(result)) {
    model <- result[[i]]
    message(names(result)[i])
    message(model$modelname)
    covariate_analysis <- length(unique(model$covariate)) > 1
    message("Covariate analysis is:", covariate_analysis)
    if (!covariate_analysis) {
      if (!grepl("model_averaging", names(result)[i])) {
        # Extract the relevant pieces of information excluding model averaging
        selected_models <- c(selected_models, model$modelname)
        CES <- c(CES, model$CES)
        CED <- c(CED, model$CED)
        CEDL <- c(CEDL, model$conf.int[1])
        CEDU <- c(CEDU, model$conf.int[2])
        AIC <- c(AIC, model$aic)
        log_likelihood <- c(log_likelihood, model$loglik)
        var <- c(var, model$MLE[1])
        a <- c(a, model$MLE[2])
        d <- c(d, model$MLE[4])
        response <- c(response, model$res.name)
      }

      # Special handling for 'model_averaging'
      if (grepl("model_averaging", names(result)[i]) && !is.null(result[[i]]$MA)) {
        message("Handling model averaging case")
        ma_info <- result[[i]]$MA
        ma_response <- result[[i]]$res.name
        ma_row <- c("Model averaging", "N/A", "N/A",
                    ma_info$conf.int.ma$BMDlower.ma,
                    ma_info$conf.int.ma$BMDupper.ma,
                    "N/A", "N/A", "N/A", "N/A", "N/A", ma_response)
        ma_rows[[length(ma_rows) + 1]] <- ma_row
        # Model weights
        model_names <- ma_info$Vmodelname
        model_weights[[length(model_weights) + 1]] <- data.frame('Selected Model' = model_names,
                                                                 'Response' = rep(ma_response, length(model_names)),
                                                                 weights = ma_info$Vweight$weight)
      }
    } else {
      # Handle covariates
      subgroups <- model$covar.txt
      message(subgroups)
      # Check if the covariate has a significant effect on the BMD
      # if not, extend the CED table to include a row for each subgroup
      ## TEST: we need to test this with 3+ subgroups. What if one subgroup is significant and the others are not?
      if (!grepl("model_averaging", names(result)[i])) {
        signif_covariate <- nrow(model$ced.table) > 1
        if (signif_covariate == FALSE) {
          message("Covariate does not have a significant effect on the BMD")
          model$ced.table <- data.frame('subgroup' = subgroups,
                              'BMDL' = rep(model$ced.table$BMDL, length(subgroups)),
                              'BMDU' = rep(model$ced.table$BMDU, length(subgroups)),
                              'BMD' = rep(model$ced.table$BMD, length(subgroups)))
        }
      }
      for (subgroup in subgroups) {
        message("Working on covariate: ", subgroup)
        if (!grepl("model_averaging", names(result)[i])) {
          message("For subgroup ", subgroup, " in ", names(result)[i])
          covariates <- c(covariates, subgroup)
          current_ced <- model$ced.table[model$ced.table$subgroup == subgroup, ]
          extra_info <- data.frame(names = model$text.par, values = model$MLE)
          selected_models <- c(selected_models, model$modelname)
          CES <- c(CES, model$CES)
          CED <- c(CED, current_ced$BMD)
          CEDL <- c(CEDL, current_ced$BMDL)
          CEDU <- c(CEDU, current_ced$BMDU)
          AIC <- c(AIC, model$aic)
          log_likelihood <- c(log_likelihood, model$loglik)
          var <- c(var, extra_info[extra_info$names == paste0("var-",subgroup),]$values)
          a <- c(a, extra_info[extra_info$names == paste0("a-",subgroup), ]$values)
          d <- c(d, extra_info[extra_info$names == "d-", ]$values)
          response <- c(response, model$res.name)
          message(length(covariates))
        }
        message(names(result)[i])
        message(result[[i]]$MA)

        if (grepl("model_averaging", names(result)[i]) && !is.null(result[[i]]$MA)) {
          message("Handling model averaging case for ", subgroup)
          ma_info <- result[[i]]$MA
          ma_response <- result[[i]]$res.name
          ma_row <- c("Model averaging", subgroup, "N/A", "N/A",
          ma_info$conf.int.ma[ma_info$conf.int.ma$subgroup == subgroup,]$BMDlower.ma,
          ma_info$conf.int.ma[ma_info$conf.int.ma$subgroup == subgroup,]$BMDupper.ma,
          "N/A", "N/A", "N/A", "N/A", "N/A", ma_response)
          ma_rows[[length(ma_rows) + 1]] <- ma_row

          # Model weights
          model_names <- ma_info$Vmodelname
          model_weights[[length(model_weights) + 1]] <- data.frame('Selected Model' = model_names,
                                                               'Response' = rep(ma_response, length(model_names)),
                                                               'Covariates' = rep(subgroup, length(model_names)),
                                                               weights = ma_info$Vweight$weight)
        }
      }
    }
  }
  # Combine the vectors into a dataframe
  result_df <- data.frame(
    'Selected Model' = selected_models,
    'CES' = CES,
    'CED' = CED,
    'CEDL' = CEDL,
    'CEDU' = CEDU,
    'AIC' = AIC,
    'Log-Likelihood' = log_likelihood,
    'Var' = var,
    'a' = a,
    'd' = d,
    'Response' = response
  )

  # if covariates is not empty, insert it into table:
  if (length(covariates) > 0) {
    result_df <- data.frame('Selected Model' = result_df[, 1],
                            'Covariates' = covariates,
                            result_df[, -1])
  }
  if (any(grepl("model_averaging", names(result)))) {
    # Add model averaging rows to the  results dataframe
    ma_frame <- setNames(as.data.frame(do.call(rbind, ma_rows)),
                         colnames(result_df))
    result_df <- rbind(result_df, ma_frame)

    # Bind models with their corresponding ma weights
    weights <- as.data.frame(do.call(rbind, model_weights))

    if (length(covariates) > 0) {
      result_df <- merge(result_df, weights, by = c("Selected.Model", "Response", "Covariates"), all.x = TRUE)
    } else {
      result_df <- merge(result_df, weights, by = c("Selected.Model", "Response"), all.x = TRUE) 
    }
  } else {
    result_df$weights <- NA
  }

  return(result_df)
}



f.plot.result <- function(proast_results_list,
                          output_path = NULL,
                          output_type = "svg",
                          prefix = NULL,
                          model_averaging = FALSE) {

  if (model_averaging == FALSE) {
    result <- proast_results_list[grepl("plot_result", names(proast_results_list))]
  } else {
    if (model_averaging == TRUE) {
      result <- proast_results_list[grepl("model_averaging", names(proast_results_list))]
    } else {
      stop("Invalid choice for model averaging. If you want to plot model averaging results, set it to TRUE, otherwise, set it to FALSE")
    }
  }

  # Output directory
  if (is.null(output_path)) {
    output_dir <- file.path(here::here(), "output")
  } else {
    output_dir <- file.path(output_path)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Output filename
  if (!is.null(prefix)) {
    filename <- file.path(output_dir, paste0("PROAST_", prefix))
  } else {
    filename <- file.path(output_dir, paste0("PROAST"))
  }
  for (i in seq_along(result)) {
    res <- result[[i]]
    if (output_type == "svg") {
      res$svg.plots <- TRUE
    }
    if (!is.null(output_type) && !output_type == "none") {
      message("Setting WAPP to TRUE")
      res$WAPP <- TRUE
      if (model_averaging == FALSE) {
        f.plot.gui(res, filename = filename, output_type = output_type)
      } else {
        if (model_averaging == TRUE) {
          f.boot.ma(res, filename = filename, output_type = output_type)
        } else {
         stop("Invalid choice for model averaging. If you want to plot model averaging results, set it to TRUE, otherwise, set it to FALSE")
        }
      }
    } else { # When output type isn't set
      if (model_averaging == FALSE) {
        f.plot.gui(res)
      } else {
        if (model_averaging == TRUE) {
          f.boot.ma(res)
        }
      }
    }
  }
}

  ### TODO
  # - add option to plot one model at a time?

  # loop over all elements of results list
  # for (i in seq_along(proast_results_list)) {
  #   message("Plotting ", names(proast_results_list)[[i]])
  #   proast_result <- proast_results_list[[i]]
  #   if (!names(proast_results_list)[[i]] == "model_averaging") {
  #     proast_result$plt.mns <- 1
  #     proast_result$fitted <- TRUE
  #     proast_result$WAPP <- TRUE
  #     f.plot.all(proast_result)
  #     dev.off()
  #   } else {
  #     if (names(proast_results_list)[[i]] == "model_averaging") {
  #       message("Plotting model averaging results")
  #       f.boot.ma(proast_result)
  #     }
  #   }
  # }
  # if (!output_type == "none") {
  #   grDevices::dev.off()
  # }

# f.plot.result(result_A[[1]])



f.remove.blanks <- function(vec) {
    if (exists("track")) 
        print("f.remove.blanks")
    LEVELS <- levels(as.factor(vec))
    LEVELS <- subset(LEVELS, LEVELS != "")
    vec <- factor(vec, level = LEVELS)
    if (exists("track")) 
        print("f.remove.blanks:   END")
    return(vec)
}



f.alert.full <- function() {
    "\nAttention: the AIC of the best model (minimum AIC) is more than two units larger than that of the full model. \nThis might indicate a problem in the data, in particular when the difference is much larger than two units (e.g. > 5). \n\nYou might check the following options:\n1. In real-life studies, not all experimental factors are completely randomized over all animals (experimental units), \ne.g. animals were housed in the same cage within a given dose group, or order of treatments were not randomized over individual animals. \nAnother option is that individual outlying animals distort the mean response of one or more treatment groups. \nThis may lead to fluctuations in the (mean) responses among treatment groups that are larger than expected from random sampling error, \nresulting in an AIC difference with the full model larger than 2 units.\n2. the data consist of subgroups not taken into account in the model (e.g. various studies, or two sexes) \n3. the data contain litter effects not taken into account \n4. the response in the top dose group deviates substantially from the fitted model (check the CI around the observed (mean) response); \n\nAssociated actions for each of these four options are:\n1. the greater scatter in (mean) responses will result in a wider BMD CI; normally, no further action is needed, \nas the BMD approach is relatively robust to such devations. You might check this by leaving out specific \ntreatment groups (one by one) and check if this has a major impact on the BMD CI.\n2. use the factor defining the subgroups as a covariate and re-analyse the data \n3. re-analyse the data with litter effects taken into account \n4. consider to leave out the top dose; it is not recommended to leave out two high dose groups.\n\n"
}



f.explain.marks <- function(ans.all) {
    color.txt <- c("black", "red", "green", "dark blue", "light blue", 
        "pink", "grey")
    color.txt <- c(color.txt, color.txt[1:3])
    mark.txt <- c("upward triangle", "cross", "diamond", "downward triangle", 
        "cross-square", "cross-plus", "diamond-plus", "circle-plus", 
        "double triangle", "square-plus")
    if (length(ans.all$gr.txt) == 1 || is.na(ans.all$displ.txt[1])) 
        gr.txt <- ans.all$gr.txt
    else gr.txt <- ans.all$displ.txt
    nr.points <- min(length(gr.txt), 10)
    group.txt.tmp <- gr.txt[1:nr.points]
    color.txt.tmp <- color.txt[1:nr.points]
    mark.txt.tmp <- mark.txt[1:nr.points]
    cat("\nThe colors in the plot relate to the following subgroups:\n")
    print(data.frame(color = color.txt.tmp, mark = mark.txt.tmp, 
        subgroup = group.txt.tmp))
}



f.plot.CI <- function(ans.all, sort = T, logscale = T, xx.lim = NA, RPFs = FALSE, 
    ref = NA, rm.NAs = TRUE, display_plots = TRUE) {
    if (exists("track")) 
        print("f.plot.CI")
    no.hill <- FALSE
    if (ans.all$do.MA) {
        CI <- ans.all$MA$conf.int.ma
        covar.txt <- CI[, 1]
        CI.hill <- matrix(c(NA, NA), nrow = 1)
        CI <- CI[, -1]
        no.hill <- TRUE
    }
    else {
        if (length(ans.all$HILL) == 0 || is.na(ans.all$HILL[1])) 
            no.hill <- T
        if (length(ans.all$EXP) == 0) 
            ans.all$EXP$conf.int <- ans.all$conf.int
        if (ans.all$dtype %in% c(4, 6)) {
            CI <- ans.all$conf.int
            CI.hill <- matrix(c(NA, NA), nrow = 1)
        }
        else {
            CI <- ans.all$EXP$conf.int
            if (!no.hill) 
                CI.hill <- ans.all$HILL$conf.int
            else CI.hill <- matrix(c(NA, NA), nrow = 1)
        }
        if (ans.all$cont) {
            if (length(CI[, 1]) == 1 && length(CI.hill[, 1]) > 
                1) 
                CI <- matrix(CI, ncol = 2, nrow = length(CI.hill[, 
                  1]), byrow = TRUE)
            if (length(CI[, 1]) > 1 && length(CI.hill[, 1]) == 
                1) 
                CI.hill <- matrix(CI.hill, ncol = 2, nrow = length(CI[, 
                  1]), byrow = TRUE)
        }
        covar.txt <- dimnames(CI)[[1]]
        if (is.null(covar.txt)) 
            covar.txt <- ans.all$fct2.txt
    }
    if (rm.NAs) {
        lst.na.1 <- is.na(CI[, 1])
        lst.na.2 <- is.na(CI[, 2])
        lst.na <- lst.na.1 | lst.na.2
        CI <- CI[!lst.na, ]
        if (!is.matrix(CI) && !is.data.frame(CI)) {
            if (length(CI[1, ]) == 2) 
                CI <- matrix(CI, ncol = 2)
            if (length(CI[1, ]) == 3) {
                covar.txt <- CI[, 1]
                CI <- as.matrix(CI[, 2:3])
            }
        }
        if (ans.all$cont && !no.hill) {
            lst.na.1 <- is.na(CI.hill[, 1])
            lst.na.2 <- is.na(CI.hill[, 2])
            lst.na <- lst.na.1 | lst.na.2
            CI.hill <- CI.hill[!lst.na, ]
            if (!is.matrix(CI.hill)) {
                if (length(CI.hill[1, ]) == 2) 
                  CI.hill <- matrix(CI.hill, ncol = 2)
                if (length(CI.hill[1, ]) == 3) 
                  CI.hill <- as.matrix(CI.hill[, 2:3])
            }
        }
        covar.txt <- covar.txt[!lst.na]
    }
    nr.lev <- length(CI[, 1])
    CES <- ans.all$CES
    if (!any(!is.na(CI))) {
        cat("\n No confidence intervals available\n")
        return()
    }
    CI.matr <- data.frame(covar = covar.txt, LB.E = CI[, 1], 
        UB.E = CI[, 2])
    if (!no.hill) {
        CI.matr$LB.H <- CI.hill[, 1]
        CI.matr$UB.H <- CI.hill[, 2]
    }
    max.x <- 1.1 * max(ans.all$x, na.rm = TRUE)
    max.CI <- max(CI, na.rm = T)
    if (is.infinite(max.CI)) 
        max.CI <- max(CI[is.finite(CI)]) * 10
    CI[, 2][is.na(CI[, 2])] <- max.CI
    if (!no.hill) {
        max.CI.hill <- max(CI.hill, na.rm = T)
        if (is.infinite(max.CI.hill)) 
            max.CI.hill <- max(CI.hill[is.finite(CI.hill)]) * 
                10
        CI.hill[, 2][is.na(CI.hill[, 2])] <- max.CI.hill
        max.CI <- max(max.CI, max.CI.hill)
    }
    lb <- CI[, 1]
    lb.pos <- lb[lb != 0]
    lb.min <- min(lb.pos, na.rm = TRUE)
    ub <- CI[, 2]
    ub.pos <- ub[ub != 0]
    ub.min <- min(ub.pos, na.rm = TRUE)
    min.CI <- min(lb.min/10, ub.min/10)
    if (!is.na(xx.lim[1])) 
        min.CI <- 10^xx.lim[1]
    CI[, 1][CI[, 1] == 0] <- min.CI
    if (!no.hill) {
        lb <- CI.hill[, 1]
        lb.pos <- lb[lb != 0]
        lb.min <- min(lb.pos, na.rm = TRUE)
        ub <- CI.hill[, 2]
        ub.pos <- ub[ub != 0]
        ub.min <- min(ub.pos, na.rm = TRUE)
        min.CI <- min(lb.min/10, ub.min/10, min.CI)
        CI.hill[, 1][CI.hill[, 1] == 0] <- min.CI
    }
    if (logscale) {
        CI <- log10(CI)
        if (!no.hill) 
            CI.hill <- log10(CI.hill)
        min.CI <- log10(min.CI)
        max.CI <- log10(max.CI)
    }
    max.y <- length(covar.txt)
    yyy.tmp <- nr.lev
    if (!no.hill) {
        max.y <- max.y * 1.5
        yyy.tmp <- nr.lev * 1.5
    }
    max.y <- max.y - 0.7
    yy.lim <- c(-0.3, max.y)
    if (is.na(xx.lim[1])) {
        CI.tmp <- CI
        CI.tmp[, 1][CI.tmp[, 1] == 0] <- NA
        if (!no.hill) 
            CI.tmp <- cbind(CI, CI.hill)
        xx.lim <- c(min(CI.tmp, na.rm = T), max.CI)
    }
    if (max.CI < max.x) {
        if (is.infinite(xx.lim[2])) 
            xx.lim[2] <- log10(max.x)
        if (is.na(xx.lim[2])) 
            xx.lim[2] <- log10(max.x)
    }
    else {
        if (is.infinite(xx.lim[2])) 
            xx.lim[2] <- log10(max.CI)
        if (is.na(xx.lim[2])) 
            xx.lim[2] <- log10(max.CI)
    }
    xx.lab <- paste("CED", round(CES[1], 3), sep = "-")
    if (logscale) 
        xx.lab <- paste("log10-", xx.lab)
    name.wapp <- paste("b", ans.all$yans, "cedCI", sep = "")
    if (display_plots) {
        f.graph.window(1, WAPP = ans.all$WAPP, title = "", name.wapp = name.wapp, 
        plotprefix = ans.all$plotprefix, svg.plots = ans.all$svg.plots)
        plot(1:max.y, 1:max.y, xlim = xx.lim, ylim = yy.lim, ylab = "", 
            yaxt = "n", xlab = xx.lab, type = "n")
        par(cex.main = 0.9)
        if (ans.all$do.MA) 
            title(main = "BMD confidence intervals based on MA")
        else if (!RPFs) {
            if (no.hill) 
                title(main = "BMD confidence intervals per subgroup")
            else title(main = "BMD confidence intervals \n(exponential and Hill, per subgroup)")
        }
        else if (RPFs) {
            if (no.hill) 
                title(main = "RPF confidence intervals per subgroup")
            else title(main = "RPF confidence intervals \n(exponential and Hill, per subgroup)")
            par(cex.main = 0.8)
            if (!is.na(ref)) 
                title(main = paste("\n\nreference: ", ref))
        }
    }
    if (sort && length(CI[, 1]) > 1) {
        ci.mn <- (CI[, 1] + CI[, 2])/2
        covar.txt <- covar.txt[order(ci.mn)]
        CI <- CI[order(ci.mn), ]
        if (!no.hill) {
            CI.hill <- CI.hill[order(ci.mn), ]
        }
    }
    l.ty <- 1
    for (ii in 1:nr.lev) {
        yyy.tmp <- yyy.tmp - 1
        if (CI[ii, 1] == min.CI) {
            l.ty <- 2
        }
        if (is.infinite(CI[ii, 2])) {
            l.ty <- 2
            CI[ii, 2] <- max.CI
        }
        if (display_plots) {
            lines(CI[ii, ], rep(yyy.tmp, 2), lty = l.ty)
            mtext(covar.txt[ii], 4, 1, adj = 0, las = 1, at = 1.01 * 
                yyy.tmp, cex = 0.8)
        }
        l.ty <- 1
        if (!no.hill) {
            yyy.tmp <- yyy.tmp - 0.2
            if (CI.hill[ii, 1] == min.CI) {
                l.ty <- 2
            }
            if (is.infinite(CI.hill[ii, 2])) {
                l.ty <- 2
                CI.hill[ii, 2] <- max.x
            }
            if (display_plots) {
                lines(CI.hill[ii, ], rep(yyy.tmp, 2), lty = l.ty)
            }
        }
        l.ty <- 1
    }
    if (exists("track")) 
        print("f.plot.CI  END")
    if (ans.all$WAPP) 
        dev.off()
    return(CI.matr)
}



