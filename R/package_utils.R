.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname, force = TRUE)
  
  print_ascii_art()

  }

#' This function prints ASCII art when the package is loaded
#' @export
print_ascii_art <- function() {
  cat("

                    .OMMo                     '0MWc                      
                     .kMMd                   ,KMN:                      
                      oWMXxooooooooooooooooooOWM0'                      
                      'OMMMWWWWWWWWWWWWWWWWWMMMNl                       
                       ,0WMNd;'''''''''',cd0WMNo.                       
                        .xNMNx;.       'l0WWW0:                         
                          ,xNMW0dc'';lkXMWN0l.                          
                            'oKMMMWNWMMWk;'.                            
                           .:xXMMN00XWMW0oc,            .               
                         .lKWWXxc'. .;oONMMNx,        ;OKo.             
                        ;0WW0l.        .;xKNMNd.   .:ONMXd.             
                       :XMMNxccccccccccccld0WMWk.  .kWXd'               
                      ;KMMMMMMMMMMMMMMMMMMMMMMMWd.  .,'                 
                     .dMM0c:::::::::::::::::::::,.     .'''''.          
                     .kMMd                            'ONNNNNk.         
                     .kMMx.                           .cddddd:.         
                      lWMN0OOOOOOOOOOOOOOOOOOOOOo.                      
                      .xWMMNK00000000000000XWMMX:  .:kd'                
                       .xWMNd'...........:d0WMK:   .kWMXd'              
                        .cKWWKo,.     .:kNMMNx'      :ONMO'             
                          .lONMN0dc:okXWWX0x,         .;l'              
                            .lKMMMMMMMMWk'                              
                           ,dKWMN0xdkXWMNkdc.                           
                         ,kNMNkc'.   .;dKWMWKl.                         
                       .lXMNk,          .lkKWWO'                        
                      .oNMMNkdxdddddddddxxxKWMM0,                       
                      :NMWNXXXXXXXXXXXXXXXXXXNMMk.                      
                     .xMMk,..................lXMX;                      
                     .OMMo                   ,KMWc                      
                     .OMWo                   '0MWc  
                     
                     
                     
 .S_sSSs     .S       S.    .S_sSSs      sSSs    sSSs    sSSs_sSSs     .S_sSSs    
.SS~YS%%b   .SS       SS.  .SS~YS%%b    d%%SP   d%%SP   d%%SP~YS%%b   .SS~YS%%b   
S%S   `S%b  S%S       S%S  S%S   `S%b  d%S'    d%S'    d%S'     `S%b  S%S   `S%b  
S%S    S%S  S%S       S%S  S%S    S%S  S%|     S%S     S%S       S%S  S%S    S%S  
S%S    S&S  S&S       S&S  S%S    d*S  S&S     S&S     S&S       S&S  S%S    d*S  
S&S    S&S  S&S       S&S  S&S   .S*S  Y&Ss    S&S_Ss  S&S       S&S  S&S   .S*S  
S&S    S&S  S&S       S&S  S&S_sdSSS   `S&&S   S&S~SP  S&S       S&S  S&S_sdSSS   
S&S    S&S  S&S       S&S  S&S~YSSY      `S*S  S&S     S&S       S&S  S&S~YSY%b   
S*S    d*S  S*b       d*S  S*S            l*S  S*b     S*b       d*S  S*S   `S%b  
S*S   .S*S  S*S.     .S*S  S*S           .S*P  S*S.    S*S.     .S*S  S*S    S%S  
S*S_sdSSS    SSSbs_sdSSS   S*S         sSS*S    SSSbs   SSSbs_sdSSSS  S*S    S&S  
SSS~YSSY      YSSP~YSSY    S*S         YSS'      YSSP    YSSP~YSSSSS  S*S    SSS  
                           SP                                         SP          
                           Y                                          Y             
")
}

