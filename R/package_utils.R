.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname, force = TRUE)
  
  print_ascii_art()

  }

#' This function prints ASCII art when the package is loaded
#' @export
print_ascii_art <- function() {
  cat("

                      .OMMo                   '0MWc                      
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

 ██████████                         █████████                     ███████████  
░░███░░░░███                       ███░░░░░███                   ░░███░░░░░███ 
 ░███   ░░███ █████ ████ ████████ ░███    ░░░   ██████   ████████ ░███    ░███ 
 ░███    ░███░░███ ░███ ░░███░░███░░█████████  ███░░███ ███░░███  ░██████████  
 ░███    ░███ ░███ ░███  ░███ ░███ ░░░░░░░░███░███████ ░███ ░███  ░███░░░░░███ 
 ░███    ███  ░███ ░███  ░███ ░███ ███    ░███░███░░░  ░███ ░███  ░███    ░███ 
 ██████████   ░░████████ ░███████ ░░█████████ ░░██████ ░░███████  █████   █████
░░░░░░░░░░     ░░░░░░░░  ░███░░░   ░░░░░░░░░   ░░░░░░   ░░░░░███ ░░░░░   ░░░░░ 
                         ░███                               ░███               
                         █████                              █████ version.version             
                        ░░░░░                              ░░░░░               
                                                               
      ")
}

