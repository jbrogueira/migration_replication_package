from config import σ

def get_kappa(sigma):
    
    table = {
        0.5: (21, 30),
        2:   (2.3,  2.6),
        #2:   (2.3,  2.55), #for higher progressivity only
        #2:   (2.05,  2.25), # for ϵ = 20
        4:   (1.58, 1.65),
        6:   (1.4, 1.45),
        8:   (1.31, 1.35),
        10:  (1.26, 1.30),
        16:  (1.19, 1.22),
        20:  (1.17, 1.2),
    }
    return table.get(sigma, (1.0, 1.0))

# assign at import
κl, κh = get_kappa(σ)