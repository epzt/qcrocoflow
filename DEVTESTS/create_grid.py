dimensions:
    xi_rho = 182 ;
    eta_rho = 122 ;
    xi_u = 181 ;
    eta_u = 122 ;
    xi_v = 182 ;
    eta_v = 121 ;
    xi_psi = 181 ;
    eta_psi = 121 ;
    s_rho = 30 ;          # Niveaux sigma pour la grille rho (ρ)
    s_w = 31 ;            # Niveaux sigma pour la grille w
    time = UNLIMITED ;    # Pour permettre la gestion temporelle

variables:
    double lon_rho(eta_rho, xi_rho) ;       # Coordonnées des points rho
        lon_rho:units = "degrees_east" ;
    double lat_rho(eta_rho, xi_rho) ;       # Coordonnées latitude pour rho
        lat_rho:units = "degrees_north" ;
    double lon_u(eta_u, xi_u) ;             # Coordonnées longitude pour u
        lon_u:units = "degrees_east" ;
    double lat_u(eta_u, xi_u) ;             # Coordonnées latitude pour u
        lat_u:units = "degrees_north" ;
    double lon_v(eta_v, xi_v) ;             # Coordonnées longitude pour v
        lon_v:units = "degrees_east" ;
    double lat_v(eta_v, xi_v) ;             # Coordonnées latitude pour v
        lat_v:units = "degrees_north" ;
    double lon_psi(eta_psi, xi_psi) ;       # Coordonnées longitude pour psi
        lon_psi:units = "degrees_east" ;
    double lat_psi(eta_psi, xi_psi) ;       # Coordonnées latitude pour psi
        lat_psi:units = "degrees_north" ;

    double h(eta_rho, xi_rho) ;             # Bathymétrie en mètres
        h:units = "meters" ;
        h:long_name = "Bathymetry" ;

    double f(eta_rho, xi_rho) ;             # Paramètre de Coriolis
        f:units = "second-1" ;
        f:long_name = "Coriolis parameter" ;

    double pm(eta_rho, xi_rho) ;            # Facteur d’échelle en x
        pm:units = "meter-1" ;
        pm:long_name = "inverse grid spacing in XI" ;

    double pn(eta_rho, xi_rho) ;            # Facteur d’échelle en y
        pn:units = "meter-1" ;
        pn:long_name = "inverse grid spacing in ETA" ;

    double angle(eta_rho, xi_rho) ;         # Angle de rotation de la grille
        angle:units = "radians" ;
        angle:long_name = "angle between XI axis and east" ;

    int mask_rho(eta_rho, xi_rho) ;         # Masque terrestre pour rho
        mask_rho:flag_values = "0, 1" ;
        mask_rho:flag_meanings = "land, water" ;

    int mask_u(eta_u, xi_u) ;               # Masque terrestre pour u
        mask_u:flag_values = "0, 1" ;
        mask_u:flag_meanings = "land, water" ;

    int mask_v(eta_v, xi_v) ;               # Masque terrestre pour v
        mask_v:flag_values = "0, 1" ;
        mask_v:flag_meanings = "land, water" ;

    int mask_psi(eta_psi, xi_psi) ;         # Masque terrestre pour psi
        mask_psi:flag_values = "0, 1" ;
        mask_psi:flag_meanings = "land, water" ;

    double s_rho(s_rho) ;                   # Niveaux sigma sur rho
        s_rho:long_name = "S-coordinate at RHO-points" ;
        s_rho:positive = "up" ;
        s_rho:standard_name = "ocean_s_coordinate" ;

    double s_w(s_w) ;                       # Niveaux sigma sur w
        s_w:long_name = "S-coordinate at W-points" ;
        s_w:positive = "up" ;
        s_w:standard_name = "ocean_s_coordinate" ;

    double hc ;                             # Profondeur critique
        hc:long_name = "S-coordinate parameter, critical depth" ;
        hc:units = "meter" ;

    double theta_s ;                        # Paramètre de surface
        theta_s:long_name = "S-coordinate surface control parameter" ;

    double theta_b ;                        # Paramètre de fond
        theta_b:long_name = "S-coordinate bottom control parameter" ;

    double Tcline ;                         # Profondeur de transition thermocline
        Tcline:long_name = "S-coordinate surface/bottom layer width" ;
        Tcline:units = "meter" ;

    # Dimensions temporelles (éventuellement utilisées dans les grilles dynamiques)
    double zeta(time, eta_rho, xi_rho) ;    # Hauteur de surface de la mer
        zeta:units = "meter" ;
        zeta:long_name = "sea surface height" ;

    double time(time) ;                     # Variable de temps
        time:units = "days since 0001-01-01 00:00:00" ;
        time:calendar = "gregorian" ;
        
        
