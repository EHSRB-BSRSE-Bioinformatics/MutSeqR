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
                     
                     
 ,ggg, ,ggg,_,ggg,                            ,gg,                       ,ggggggggggg,   
dP  Y8dP  Y88P  Y8b                I8        i8  8i                     dP   88      Y8, 
Yb, `88'  `88'  `88                I8        `8,,8'                     Yb,  88      `8b 
 `   88    88    88             88888888      `88'                       `   88      ,8P 
     88    88    88                I8         dP 8,                          88aaaad8P   
     88    88    88  gg      gg    I8        dP' `8a   ,ggg,     ,gggg,gg    88    Yb,   
     88    88    88  I8      II    I8       dP'   `Yb i8   8i   dP    Y8I    88      8b  
     88    88    88  I8,    ,8I   ,I8,  _ ,dP'     I8 I8, ,8I  i8'    ,8I    88      `8i 
     88    88    Y8,,d8b,  ,d8b, ,d88b,  888,,____,dP `YbadP' ,d8,   ,d8b    88       Yb,
     88    88    `Y88P' Y88P `Y888P  Y88a8P Y88888P  888P Y888P Y8888P 88d   88        Y8
                                                                       I8P               
                                                                       I8'               
                                                                       I8                
                                                                       I8                
                                                                       I8                                    

")
}

