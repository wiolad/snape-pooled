let rho epsa epsr = epsr *. (1. -. 2. *. epsa) /. (1. -. epsa -. epsr);;
let alpha  epsa epsr = epsa *. (1. -. 2. *. epsr) /. ( 1. -. epsa -. epsr);;

let ptilde epsa epsr p =
    let rho = rho epsa epsr and
    alpha = alpha epsa epsr 
    and q = 1. -. p in
    (1. -. rho) *. p +. alpha *. q 
;;

let epsprime eps = ( 1. -. eps ) /. eps 
;;

let x epsa epsr p =
    let q = 1. -. p and
    ea = epsprime epsa
    and er = epsprime epsr
    in ea /. (1. -. ea*.er) *. (q /. p -. er);;

let y epsa epsr p =
let q = 1. -. p and
    ea = epsprime epsa
    and er = epsprime epsr
    in er /. (1. -. ea*.er) *. (p /. q -. ea);;

let ptilde2 epsa epsr p = 
    let q = 1. -. p and
    x = x epsa epsr p
    and y = y epsa epsr p
    in x *. p +. (1. -. y) *. q
;;

let rec range start stop acc =
    if (start > stop) then
        List.rev acc
    else
        range (start +1) stop (start::acc)
;;

let p = List.map (fun e -> float_of_int e /. 10.0) (range 1 9 []) 
    and errors = [ 1e-6; 1e-5; 1e-3 ; 1e-4; 1e-2; 1e-1 ] in
    let epsar = List.flatten ( List.map ( fun f -> List.map 
                    (fun e -> (f , e ) ) errors ) errors ) 
    in
    Printf.fprintf stderr "1:p 2:ea 3:er 4:ptilde 5:ptilde2 6:x 7:y 8:diff 9:abs_diff\n"; 
    List.iter (fun t ->
    List.iter (fun e  ->
    let epsa = fst t and epsr =snd t in
    let ptilde= ptilde epsa epsr e 
    and ptilde2 = ptilde2 epsa epsr e in
    let diff = (ptilde -. ptilde2) in
    let abs_diff =  ( abs_float diff ) in
    Printf.fprintf stdout 
        "%g %g %g %g %g %g %g %g %g\n" 
                e epsa epsr 
                ptilde ptilde2 
                (x epsa epsr e ) (y epsa epsr e)
                diff abs_diff    
                ) p 
                ) epsar 
;;
