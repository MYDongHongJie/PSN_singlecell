
	$('#reportTable').bootstrapTable({
		height: 300,//定义整个表格的高度
		pagination: true,//显示分页
		pageSize: 10,//设置每页多少条数据，默认是10条
		pageNumber: 1,//默认第一页
		showColumns: true,//显示隐藏选择几个数据类型的按钮
		showExport: true,//显示导出文件的按钮
		exportDataType: 'basic',//basic(默认), all(导出全部)
		exportTypes: ['csv', 'xls', 'doc', 'txt'],//定义下载文档格式
		search: true,//显示搜索框
		tableTitle: '表2 下机数据统计',//定义表格的标题
		optionalParams: {
			fileName: '下机数据统计',//自定义文件名，默认是tableExport
		},//自定义title内容(鼠标移到表格的thead上显示的内容)
		titleContent: ['Sample_Name', 'Read_Num', 'Total_base', 'N_rate', 'GC_Content', 'Q20_rate', 'Q30_rate', 'Average_Read_Num', 'N_num', 'GC_Num', 'Q20_Num', 'Q30_Num', 'Read_Num', 'Total_base', 'N_rate', 'GC_Content', 'Q20_rate', 'Q30_rate', 'Average_Read_Num', 'N_num', 'GC_Num', 'Q20_Num', 'Q30_Num'],
		columns: tableColumns,//表格的thead数据
		data: tableData,//表格的tbody数据
		ellipsis: false,//超出隐藏thead的超出内容
		formatLoadingMessage: function () {//加载中触发这个方法
			return "请稍等，正在加载中...";
		},
		onSearch: function (text) {
		},
		onPageChange: function (size, number) {//分页改变时触发这个方法

		},
		formatNoMatches: function() {//搜索找不到内容时触发这个方法
			return '暂无相关内容！';
		}
	});
	//窗口变化表格随着改变						
	$(window).resize(function () {
		$('#reportTable').bootstrapTable('resetView');
	});

	var tableCont = document.querySelector('.fixed-table-body');
	// var tableLeft = document.querySelector('.fix-table-sidebar tbody');
	function scrollHandle(e) {
		var scrollTop = this.scrollTop;
		this.querySelector('thead').style.transform = 'translateY(' + scrollTop + 'px)';
		// tableLeft.style.transform = 'translateY(' + -scrollTop + 'px)';
	}
	
	// tableCont.addEventListener('scroll', scrollHandle);

	function slideNav() {
		var top = $(document).scrollTop();
		if (top >= 62) {
			$('#slide-nav').removeClass('slideDown').addClass('slideUp');
			$('#sidebar, .content-id1, .content-id2, .content-id3, .content-id4, .catalog-id1, .catalog-id2, .catalog-id3, .catalog-id4').removeClass('positiontop').addClass('positionbottom');
		} else {
			$('#slide-nav').removeClass('slideUp').addClass('slideDown');
			$('#sidebar, .content-id1, .content-id2, .content-id3, .content-id4, .catalog-id1, .catalog-id2, .catalog-id3, .catalog-id4').removeClass('positionbottom').addClass('positiontop');
		}

	}

	slideNav();
	$(document).scroll(function () {
		slideNav();
	})
	
	var lBar = '';
	//样品数组
	var carouselArr = ['样品1', '样品2', '样品3', '样品4', '样品5', '样品6', '样品7', '样品8', '样品9', '样品10', '样品11', '样品12','样品13', '样品14', '样品15', '样品16'];
	//样品循环
	carouselArr.forEach(function (val, index) {
		lBar += 
				'<li class="'+ (index == 0 ? 'cl-list-active': '') +'">'+ val +'</li>';
	})
	$('.cl-list').html(lBar);

	//点击切换样品
	$('.cl-list li').click(function () {
		var index = $(this).index();
		$('.cl-list li').removeClass('cl-list-active');
		$(this).addClass('cl-list-active');
		$('.cr-img, .carousel-img-zoom img').attr('src', './result/img/'+ (index + 1) +'.png');
		$('.cr-download-ul-a').attr('href', './result/img/'+ (index + 1) +'.png');
	})

	//模糊搜索匹配
	var searchVal = function (str) {
		var arr = [];
		str = str.replace(/(^\s*)|(\s*$)/g, '');
		carouselArr.forEach(function (val, index) {
			if (val.indexOf(str) > -1) {
				arr.push(val);
			}
		})
		if (str.length == 0) {
			$('.cl-search-result').hide();
		} else if (arr.length == 0) {
			$('.cl-search-result').show().html('<li>暂无相关内容</li>');
		} else {
			var strLi = '';
			arr.forEach(function (val) {
				strLi += '<li>'+ val +'</li>';
			})
			$('.cl-search-result').show().html(strLi);
		}
	};
	//搜索结果
	$('.cl-search-in').on('keyup', function () {
		setTimeout(function () {
			var val = $('.cl-search-in').val();
			searchVal(val);
		}, 500)
		
	})

	//选择样品
	$('.cl-search-result').on('click', 'li', function () {
		var txt = $(this).html();
		if (txt == '暂无相关内容')return;
		$('.cl-search-in').val(txt);
		$('.cl-search-result').hide();
		var index = carouselArr.indexOf(txt);
		$('.cl-list li').removeClass('cl-list-active');
		$('.cl-list li:eq('+ index +')').addClass('cl-list-active');
		$('.cr-img, .carousel-img-zoom img').attr('src', './result/img/'+ (index + 1) +'.png');
		$('.cr-download-ul-a').attr('href', './result/img/'+ (index + 1) +'.png');
	})


	//点击搜索按
	$('.cl-search-icon').click(function () {
		var index = carouselArr.indexOf($('.cl-search-in').val());
		if (index == -1) return;
		$('.cl-search-result').hide();
		$('.cl-list li').removeClass('cl-list-active');
		$('.cl-list li:eq('+ index +')').addClass('cl-list-active');
		$('.cr-img, .carousel-img-zoom img').attr('src', './result/img/'+ (index + 1) +'.png');
		$('.cr-download-ul-a').attr('href', './result/img/'+ (index + 1) +'.png');
		
	});


	//显示隐藏下载图片cc-download-active
	$('.cr-download-btn').click(function () {
		var active = $('.cr-download-btn').hasClass('cr-download-active');
		if (active) {
			$('.cr-download-ul').hide();
			$('.cr-download-btn').removeClass('cr-download-active');
		} else {
			$('.cr-download-ul').show();
			$('.cr-download-btn').addClass('cr-download-active');
		}
	})

	//图片放大
	$('.cr-img').click(function () {
		$('.carousel-img-zoom').show();
	})
	//图片缩放
	$('.carousel-img-zoom img').click(function () {
		$('.carousel-img-zoom').hide();
	})

	//鼠标移上背景变色
	$('.cl-list li').mouseenter(function () {
		if (!$(this).hasClass('cl-list-active')) {
			$(this).addClass('cl-list-hover');
		}
	}).mouseleave(function () {
		$(this).removeClass('cl-list-hover');
	})

	//图片和左侧栏随窗口变化
    
	function resizWidth() {
		var winWidth = $('body').width();
		if (winWidth <= 1200)
			winWidth = 1200;
		var leftHeight = parseInt((462 / 1423) * winWidth);
		var imgWidth = parseInt((613 / 1423) * winWidth);
		var maxHeight = parseInt((392 / 1423) * winWidth);
		$('.cl-list').css('max-height', maxHeight + 'px');
		$('.carousel-left').css('height', leftHeight + 'px');
		$('.cr-img').css('width', imgWidth + 'px');
		$('.carousel-right').css('width', (imgWidth + 1) + 'px');
	}
	resizWidth();

	//判断表格的x轴是否有滚动条
	function tableScroll() {
		var tableEle = $('#reportTable').parent('.fixed-table-body');
		if (tableEle.height() - 1 >= $('#reportTable').height()) {
			$('.fix-table-sidebar').css({
				'overflow': 'auto',
				'height': 'auto'
			});
		} else {
			if (tableEle.css('overflow-x') == 'hidden' || tableEle.width() > $('#reportTable').width()) {
				$('.fix-table-sidebar').css({
					'overflow': 'hidden',
					'height': (tableEle.height() - 2) + 'px'
				});
			} else {
				$('.fix-table-sidebar').css({
					'overflow': 'hidden',
					'height': (tableEle.height() - 19) + 'px'
				});
			}
		}
	}

	// tableScroll();

	$(window).resize(function () {
		resizWidth();
		tableScroll();
	})

	//记录每个滚动条离开时的位置
	var scrollObj = {
		proTop: 0,
		baseTop: 0,
		rnaTop: 0,
		appTop: 0
	};

	//返回上次的位置
	var recordName = 'proTop';
	$('.justify-content-center').on('click', 'li', function () {
		scrollObj[recordName] = $(window).scrollTop();
		if ($(this).hasClass('justify-item1')) {
			recordName = 'proTop';
		} else if ($(this).hasClass('justify-item2')) {
			recordName = 'baseTop';
		} else if ($(this).hasClass('justify-item3')) {
			recordName = 'rnaTop';
		} else if ($(this).hasClass('justify-item4')) {
			recordName = 'appTop';
		}
		$(window).scrollTop(scrollObj[recordName]);
	})


	//根据滚动位置给左侧栏的标题添加颜色
	function resizeTop(arg, name, typeName, top) {
		$(name + ' a').removeClass('active-color');
		for (var i = 0; i < arg; i++) {
			if (i == arg - 1) {
				$(name + ' a:eq('+ i +')').addClass('active-color');
			} else {
				if (top < $(typeName + (i + 2)).offset().top - 100) {
					$(name + ' a:eq('+ i +')').addClass('active-color');
					break;
				}
			}
		}
	}


	$(window).scroll(function () {
		var top = $(window).scrollTop();
		moveUp(top);
		if ($('.justify-item1 a').hasClass('active')) {
			resizeTop(4, '.nav-pills-one', '.pro-item', top);
		} else if ($('.justify-item2 a').hasClass('active')) {
			resizeTop(11, '.nav-pills-two', '.base-item', top);
		} else if ($('.justify-item3 a').hasClass('active')) {
			resizeTop(20, '.nav-pills-three', '.rna-item', top);
		} else if ($('.justify-item4 a').hasClass('active')) {
			resizeTop(20, '.nav-pills-four', '.app-item', top);
		}
	})

	//置顶滚动按钮显示隐藏
	function moveUp(top) {
		if (top) {
			if (top > 100) {
				$('#floatbar_goTop').show('slow');
				$('.tel-box').css('margin-top', '17px');
			} else {
				$('#floatbar_goTop').hide('slow');
				$('.tel-box').css('margin-top', '-17px');
			}
		}
	}

	//展开收起悬浮框
	$('.packup').click(function () {
		if ($('.packup').hasClass('packup-active')) {
			$('.packup').removeClass('packup-active');
			$('#floatbar_goTop p, .floatbar-wrap').hide('slow');
			$('.packup i').css('transform', 'rotate(360deg)');
		} else {
			$('.packup').addClass('packup-active');
			$('#floatbar_goTop p, .floatbar-wrap').show('slow');
			$('.packup i').css('transform', 'rotate(-45deg)');
		}
	})


	//联系方式移入移出
	var timer = null;
	$('.float_btn').mouseenter(function () {
		clearTimeout(timer);
		$('.tel-box').show();
	})
	.mouseleave(function () {
		timer = setTimeout(function () {
			$('.tel-box').hide();
		}, 500)
	})
	//电话移入移出
	$('.tel-box').mouseenter(function () {
		clearTimeout(timer);
		$('.tel-box').show();
	})
	.mouseleave(function () {
		$('.tel-box').hide();
	})



	


	

	

