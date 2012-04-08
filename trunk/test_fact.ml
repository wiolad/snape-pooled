open Random

let memoize f =
    let t = Hashtbl.create 1000  in 
    fun n ->
		try  Hashtbl.find t n 
    with Not_found ->    
    (* Printf.fprintf stdout  "new2%!"; *) 
        let res = f n  in
        (Hashtbl.add t n res;
        res) 
;; 
let gammaln (z : float) =
    let cof =[|76.18009172947146;-86.50532032941677;
              24.01409824083091;-1.231739572450155;
              0.1208650973866179e-2;-0.5395239384953e-5|]
    in let y = ref z in
    let x = z in
    let tmp=x +.5.5 in
    let tmp = tmp -. (x +. 0.5)*.log(tmp) in
    let ser= ref 1.000000000190015 in
    for j= 0 to 5  do 
        y:=!y +. 1.0;
        ser := !ser +. cof.(j) /. !y 
    done;
-. tmp +. log(2.5066282746310005*. !ser /. x)
;;
let factln  =
    memoize (fun (n:int) -> gammaln (float_of_int n +. 1.0)) 
;;    
let _ = for i = 0 to 1000000 do
	let r = Random.int 100 in
	Printf.fprintf stdout "%d:%f\n" r (exp (factln r))
	done
;;
