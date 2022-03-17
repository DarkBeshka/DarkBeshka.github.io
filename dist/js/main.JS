function moveGame(e){
	var object1=document.getElementById("minigame");
	var cs=window.getComputedStyle(minigame);
	var left=parseInt(cs.marginLeft);
	var top=parseInt(cs.marginTop);

	switch(e.key){

		case "ArrowLeft":
			if (left>0)
				object1.style.marginLeft=left-10+"px";
			break;
		case "ArrowUp":
			if(top>0)
				object1.style.marginTop=top-10+"px";
			break;
		case "ArrowRight":
			if (left<document.documentElement.clientWidth - 100)
				object1.style.marginLeft=left+10+"px";
			break;
		case "ArrowDown":
			if (top<document.documentElement.clientHeight - 100)
				object1.style.marginTop=top+10+"px";
			break;
	}
}

addEventListener("keydown", moveGame);