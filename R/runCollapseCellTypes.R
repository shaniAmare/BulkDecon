
runCollapseCellTypes <- function(object, matching = NULL){
              
              if(is.null(matching)){
                  stop("matching must be set")
              }
              
              fit <- list(beta=t(Biobase::pData(object)$beta),
                          p=t(Biobase::pData(object)$p),
                          t=t(Biobase::pData(object)$t),
                          se=t(Biobase::pData(object)$se),
                          prop_of_all=t(Biobase::pData(object)$prop_of_all),
                          prop_of_nontumor=t(Biobase::pData(object)$prop_of_nontumor),
                          sigma=Biobase::pData(object)$sigmas)
              
              temp <- collapseCellTypes(fit = fit, 
                                        matching = matching)
              
              Biobase::pData(object)$beta <- t(temp$beta)
              Biobase::pData(object)$p <- t(temp$p)
              Biobase::pData(object)$t <- t(temp$t)
              Biobase::pData(object)$se <- t(temp$se)
              Biobase::pData(object)$prop_of_all <- t(temp$prop_of_all)
              Biobase::pData(object)$prop_of_nontumor <- t(temp$prop_of_nontumor)
              Biobase::pData(object)$sigmas <- temp$sigma
              
              return(object)
          }
