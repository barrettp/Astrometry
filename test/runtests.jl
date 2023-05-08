using Astrometry
using Test

@testset "Astrometry.jl" begin
    
    #   Test IAU 1980 Precession-Nutation Model
    @test all(abs.(iau_1980_nutation(2453736.5) .- (-9.643658353226563966e-6, 4.060051006879713322e-5)) .<= 1e-13)
    @test     abs(iau_1980_obliquity(2454388.5) - 0.4090751347643816218) .<= 1e-13
    @test all(abs.(iau_1976_precession(2433282.5, 2451544.5) .- 
                   (0.5588961642000161243e-2, 0.5589922365870680624e-2, 0.4858945471687296760e-2)) .<= 1e-13)
    #   nutation matrix
    rtruen80 = [0.9999999999534999268 0.8847935789636432161e-5 0.3835906502164019142e-5;
                -0.8847780042583435924e-5 0.9999999991366569963 -0.4060052702727130809e-4;
                -0.3836265729708478796e-5 0.4060049308612638555e-4 0.9999999991684415129]
    date = 2453736.5
    ϵA = iau_1980_obliquity(date)
    ψ, ϵ = iau_1980_nutation(date)
    @test all(abs.(Astrometry.Rx(-(ϵA + ϵ))*Astrometry.Rz(-ψ)*Astrometry.Rx(ϵA) .- rtruen80) .<= 1e-13)
    #   precession matrix
    rtruep76 = [0.9999995504328350733 0.8696632209480960785e-3 0.3779153474959888345e-3;
                -0.8696632209485112192e-3 0.9999996218428560614 -0.1643284776111886407e-6;
                -0.3779153474950335077e-3 -0.1643306746147366896e-6 0.9999999285899790119]
    ψ, ω, χ = iau_1976_precession(JD2000, 2450124.4999)
    @test all(abs.(Astrometry.Rz(-ω)*Astrometry.Ry(χ)*Astrometry.Rz(-ψ) .- rtruep76) .<= 1e-13)
    #   precession-nutation matrix
    rtruepn80 = [0.9999995831934611169     0.8373654045728124011e-3 0.3639121916933106191e-3;
                 -0.8373804896118301316e-3  0.9999996485439674092    0.4130202510421549752e-4;
                 -0.3638774789072144473e-3 -0.4160674085851722359e-4 0.9999999329310274805]
    @test all(abs.(precession_nutation(2450124.4999; model=:iau1980) .- rtruepn80) .<= 1e-13)

    #   Test IAU 2000 Equinox Model
    date = 2453736.5
    @test     abs(iau_2000_equinox_complement(date) - 0.2046085004885125264e-8) <= 1e-13

    #   Earth position and velocity (DE405 ephemeris)
    truehpos = [-0.7757238809297706813, 0.5598052241363340596, 0.2426998466481686993]
    truehvel = [-0.1091891824147313846e-1, -0.1247187268440845008e-1, -0.5407569418065039061e-2]
    date = 2453412.02501161
    @test all(abs.(iau_2000_position(date; frame=:heliocenter)[1] .- truehpos) .<= 1e-11)
    @test all(abs.(iau_2000_position(date; frame=:heliocenter)[2] .- truehvel) .<= 1e-13)

    truebpos = [-0.7714104440491111971, 0.5598412061824171323, 0.2425996277722452400]
    truebvel = [-0.1091874268116823295e-1, -0.1246525461732861538e-1, -0.5404773180966231279e-2]
    @test all(abs.(iau_2000_position(date; frame=:barycenter)[1] .- truebpos) .<= 1e-11)
    @test all(abs.(iau_2000_position(date; frame=:barycenter)[2] .- truebvel) .<= 1e-13)

    #   Test IAU 2000A and 2000B Precession-Nutation Model
    @test all(abs.(iau_2000a_nutation(2453736.5) .- (-0.9630909107115518431e-5, 0.4063239174001678710e-4)) .<= 1e-13)
    @test all(abs.(iau_2000b_nutation(2453736.5) .- (-0.9632552291148362783e-5, 0.4063197106621159367e-4)) .<= 1e-13)
    @test all(abs.(iau_2000_precession(2453736.5) .- (0.0014656153365596443, 0.4090927977670502, 3.028075726363971e-6)) .<= 1e-13)
    #   nutation matrix
    rtrue = [0.9999999999536227949 0.8836238544090873336e-5 0.3830835237722400669e-5;
             -0.8836082880798569274e-5 0.9999999991354655028 -0.4063240865362499850e-4;
             -0.3831194272065995866e-5 0.4063237480216291775e-4 0.9999999991671660338]
    date = 2453736.5
    Δt = (date - JD2000)/(100*Astrometry.DAYPERYEAR)
    ϵA = iau_1980_obliquity(date)+ deg2rad(1/3600)*Astrometry.ϵ_corr_2000*Δt
    ψn, ϵn = iau_2000a_nutation(date)
    @test all(abs.(Astrometry.Rx(-(ϵA + ϵn))*Astrometry.Rz(-ψn)*Astrometry.Rx(ϵA) .- rtrue) .<= 1e-13)
    #   nutation matrix
    rtrue = [0.9999999999536069682 0.8837746144871248011e-5 0.3831488838252202945e-5;
             -0.8837590456632304720e-5 0.9999999991354692733 -0.4063198798559591654e-4;
             -0.3831847930134941271e-5 0.4063195412258168380e-4 0.9999999991671806225]
    date = 2453736.5
    Δt = (date - JD2000)/(100*Astrometry.DAYPERYEAR)
    ϵA = iau_1980_obliquity(date) + deg2rad(1/3600)*Astrometry.ϵ_corr_2000*Δt
    ψn, ϵn = iau_2000b_nutation(date)
    @test all(abs.(Astrometry.Rx(-(ϵA + ϵn))*Astrometry.Rz(-ψn)*Astrometry.Rx(ϵA) .- rtrue) .<= 1e-13)
    #   precession matrix
    rtrueb  = [0.9999999999999942498 -0.7078279744199196626e-7 0.8056217146976134152e-7;
               0.7078279477857337206e-7 0.9999999999999969484 0.3306041454222136517e-7;
               -0.8056217380986972157e-7 -0.3306040883980552500e-7 0.9999999999999962084]
    rab, ψb, ϵb, ϵ0 = deg2rad(1/3600).*(Astrometry.icrs_ra_2000, Astrometry.ψ_bias_2000,
                                        Astrometry.ϵ_bias_2000, Astrometry.ϵ0_2000)
    rcalcb = Astrometry.Rx(-ϵb)*Astrometry.Ry(ψb*sin(ϵ0))*Astrometry.Rz(rab)
    @test all(abs.(rcalcb .- rtrueb) .<= 1e-13)
    #
    rtruep  = [0.9999989300532289018 -0.1341647226791824349e-2 -0.5829880927190296547e-3;
               0.1341647231069759008e-2 0.9999990999908750433 -0.3837444441583715468e-6;
               0.5829880828740957684e-3 -0.3984203267708834759e-6 0.9999998300623538046]
    ψ, ω, χ = iau_2000_precession(date)
    rcalcp = Astrometry.Rz(χ)*Astrometry.Rx(-ω)*Astrometry.Rz(-ψ)*Astrometry.Rx(ϵ0)
    @test all(abs.(rcalcp .- rtruep) .<= 1e-13)
    #
    rtruebp = [0.9999989300052243993 -0.1341717990239703727e-2 -0.5829075749891684053e-3;
               0.1341718013831739992e-2 0.9999990998959191343 -0.3505759733565421170e-6;
               0.5829075206857717883e-3 -0.4315219955198608970e-6 0.9999998301093036269]
    rcalcbp = (Astrometry.Rz(χ)*Astrometry.Rx(-ω)*Astrometry.Rz(-ψ)*Astrometry.Rx(ϵ0) *
              Astrometry.Rx(-ϵb)*Astrometry.Ry(ψb*sin(ϵ0))*Astrometry.Rz(rab))
    @test all(abs.(rcalcbp .- rtruebp) .<= 1e-13)
    #
    rtruen  = [0.9999999999536069682 0.8837746144872140812e-5 0.3831488838252590008e-5;
               -0.8837590456633197506e-5 0.9999999991354692733 -0.4063198798559573702e-4;
               -0.3831847930135328368e-5 0.4063195412258150427e-4 0.9999999991671806225]
    rcalcn  = Astrometry.Rx(-(ϵA + ϵn))*Astrometry.Rz(-ψn)*Astrometry.Rx(ϵA)
    @test all(abs.(rcalcn .- rtruen) .<= 1e-13)
    #
    rtruepn = [0.9999989440499982806 -0.1332880253640848301e-2 -0.5790760898731087295e-3;
               0.1332856746979948745e-2 0.9999991109064768883 -0.4097740555723063806e-4;
               0.5791301929950205000e-3 0.4020553681373702931e-4 0.9999998314958529887]
    rcalcpn = (Astrometry.Rx(-(ϵA + ϵn))*Astrometry.Rz(-ψn)*Astrometry.Rx(ϵA) *
               Astrometry.Rz(χ)*Astrometry.Rx(-ω)*Astrometry.Rz(-ψ)*Astrometry.Rx(ϵ0) *
               Astrometry.Rx(-ϵb)*Astrometry.Ry(ψb*sin(ϵ0))*Astrometry.Rz(rab))
    @test all(abs.(rcalcpn .- rtruepn) .<= 1e-13)    
    #   precession-nutation matrix
    rtrue00a = [0.9999989440476103435 -0.1332881761240011763e-2 -0.5790767434730085751e-3;
                0.1332858254308954658e-2 0.9999991109044505577 -0.4097782710396580452e-4;
                0.5791308472168152904e-3 0.4020595661591500259e-4 0.9999998314954572304]
    date = 2453736.5
    @test all(abs.(precession_nutation(date; model=:iau2000a) .- rtrue00a .<= 1e-13))
    #   precession-nutation matrix
    rtrue00b = [0.9999989440499982806 -0.1332880253640849194e-2 -0.5790760898731091166e-3;
                0.1332856746979949638e-2 0.9999991109064768883 -0.4097740555723081811e-4;
                0.5791301929950208873e-3 0.4020553681373720832e-4 0.9999998314958529887]
    date = 2453736.5
    @test all(abs.(precession_nutation(date; model=:iau2000b) .- rtrue00b) .<= 1e-13)

    #   Test IAU 2006 Equinox Model
    date = 2453736.5
    coord = (0.5791308486706011000e-3, 0.4020579816732961219e-4)
    @test abs(iau_2006_cio_locator(date, coord) - -0.1220032213076463117e-7) .<= 1e-13
        
    
    #   Test IAU 2006 Precession-Nutation Model
    date = 2453736.5
    @test all(abs.(iau_2006_cip_xy(date) .- (0.5791308486706010975e-3, 0.4020579816732958141e-4)) .<= 1e-13)
    
    truen06 = (-0.9630912025820308797e-5, 0.4063238496887249798e-4)
    date = 2453736.5
    @test all(abs.(iau_2006a_nutation(date) .- truen06) .<= 1e-13)
    #
    date = 2454388.5
    @test     abs(iau_2006_obliquity(date) - 0.4090749229387258204) .<= 1e-13
    #
    truep06 = (-0.2243387670997995690e-5, 0.4091014602391312808, -0.9501954178013031895e-3, 0.4091014316587367491)
    date = 2450124.4999
    @test all(abs.(iau_2006_precession(date) .- truep06) .<= 1e-13)
    #   precession nutation matrix
    truepn06 = [0.9999995832794205484 0.8372382772630962111e-3 0.3639684771140623099e-3;
                -0.8372533744743683605e-3 0.9999996486492861646 0.4132905944611019498e-4;
                -0.3639337469629464969e-3 -0.4163377605910663999e-4 0.9999999329094260057]
    date = 2450124.4999
    @test all(abs.(precession_nutation(date; model=:iau2006a) .- truepn06) .<= 1e-13)

    #   Barycentric Dynamical Time minus Terrestrial Time
    @test     abs(iau_2006_tdb_tt(2448939.623, 0.76543, 5.0123, 5525.242, 3190.0) - -0.1280368005936998991e-2) <= 1e-13

    #   Target direction unit vector in BCRS
    truepos = (0.2328137623960308438, 0.6651097085397855328, 0.7095257765896359837)
    position, pmotion, parallax, rvelocity = (1.234, 0.789), (1e-5, -2e-5), 1e-2, 10.0
    pmt, observer = 8.75, (0.9, 0.4, 0.1)
    @test all(abs.(proper_motion(position, pmotion, parallax, rvelocity, pmt, observer) .- truepos) .<= 1e-13)
end


