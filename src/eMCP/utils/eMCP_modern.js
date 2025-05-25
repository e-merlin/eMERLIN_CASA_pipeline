// Simple JavaScript enhancements for eMCP weblog
// Provides basic interactivity without requiring additional libraries

// Wait for DOM to be fully loaded
document.addEventListener('DOMContentLoaded', function() {
    // Mobile navigation toggle
    const navToggle = document.querySelector('.nav-toggle');
    const navList = document.querySelector('.nav-list');
    
    if (navToggle) {
        navToggle.addEventListener('click', function() {
            navList.classList.toggle('show');
        });
    }
    
    // Make tables sortable
    initSortableTables();
    
    // Initialize image zoom functionality
    initImageZoom();
    
    // Initialize any collapsible sections
    initCollapsible();
});

// Simple table sorting functionality
function initSortableTables() {
    const tables = document.querySelectorAll('.sortable');
    
    tables.forEach(table => {
        const headers = table.querySelectorAll('th');
        
        headers.forEach((header, index) => {
            if (!header.classList.contains('no-sort')) {
                header.style.cursor = 'pointer';
                header.addEventListener('click', function() {
                    sortTable(table, index);
                });
                
                // Add sort indicator
                const span = document.createElement('span');
                span.className = 'sort-indicator';
                span.textContent = ' ↕';
                header.appendChild(span);
            }
        });
    });
}

// Sort table by column
function sortTable(table, column) {
    const rows = Array.from(table.querySelectorAll('tbody tr'));
    const headers = table.querySelectorAll('th');
    const currentHeader = headers[column];
    
    // Determine sort direction
    const isAscending = currentHeader.classList.contains('sort-asc');
    
    // Reset all headers
    headers.forEach(header => {
        header.classList.remove('sort-asc', 'sort-desc');
        const indicator = header.querySelector('.sort-indicator');
        if (indicator) indicator.textContent = ' ↕';
    });
    
    // Set new sort direction
    if (isAscending) {
        currentHeader.classList.add('sort-desc');
        currentHeader.querySelector('.sort-indicator').textContent = ' ↓';
    } else {
        currentHeader.classList.add('sort-asc');
        currentHeader.querySelector('.sort-indicator').textContent = ' ↑';
    }
    
    // Sort rows
    rows.sort((a, b) => {
        let aValue = a.cells[column].textContent.trim();
        let bValue = b.cells[column].textContent.trim();
        
        // Try to convert to number if possible
        const aNum = parseFloat(aValue);
        const bNum = parseFloat(bValue);
        
        if (!isNaN(aNum) && !isNaN(bNum)) {
            return isAscending ? bNum - aNum : aNum - bNum;
        } else {
            return isAscending ? 
                bValue.localeCompare(aValue) : 
                aValue.localeCompare(bValue);
        }
    });
    
    // Reorder rows in the table
    const tbody = table.querySelector('tbody');
    rows.forEach(row => tbody.appendChild(row));
}

// Image zoom functionality
function initImageZoom() {
    const images = document.querySelectorAll('.zoomable-img');
    
    images.forEach(img => {
        img.style.cursor = 'zoom-in';
        
        img.addEventListener('click', function() {
            const modal = document.createElement('div');
            modal.className = 'image-modal';
            modal.style.position = 'fixed';
            modal.style.top = '0';
            modal.style.left = '0';
            modal.style.width = '100%';
            modal.style.height = '100%';
            modal.style.backgroundColor = 'rgba(0,0,0,0.8)';
            modal.style.zIndex = '1000';
            modal.style.display = 'flex';
            modal.style.alignItems = 'center';
            modal.style.justifyContent = 'center';
            
            const modalImg = document.createElement('img');
            modalImg.src = this.src;
            modalImg.style.maxWidth = '90%';
            modalImg.style.maxHeight = '90%';
            modalImg.style.objectFit = 'contain';
            
            modal.appendChild(modalImg);
            document.body.appendChild(modal);
            
            modal.addEventListener('click', function() {
                document.body.removeChild(modal);
            });
        });
    });
}

// Collapsible sections
function initCollapsible() {
    const collapsibles = document.querySelectorAll('.collapsible-header');
    
    collapsibles.forEach(header => {
        header.style.cursor = 'pointer';
        
        // Add indicator
        const indicator = document.createElement('span');
        indicator.className = 'collapse-indicator';
        indicator.textContent = ' ▼';
        header.appendChild(indicator);
        
        const content = header.nextElementSibling;
        content.classList.add('collapsible-content');
        
        header.addEventListener('click', function() {
            const isCollapsed = content.style.display === 'none';
            
            content.style.display = isCollapsed ? 'block' : 'none';
            indicator.textContent = isCollapsed ? ' ▼' : ' ►';
        });
    });
}

// Simple image comparison slider
function initComparisonSlider() {
    const sliders = document.querySelectorAll('.img-comp-container');
    
    sliders.forEach(container => {
        const width = container.offsetWidth;
        const slider = container.querySelector('.img-comp-slider');
        const beforeImg = container.querySelector('.img-comp-before');
        
        // Set initial position
        slider.style.left = (width / 2) + 'px';
        beforeImg.style.width = (width / 2) + 'px';
        
        // Slider functionality
        slider.addEventListener('mousedown', startSlide);
        
        function startSlide(e) {
            e.preventDefault();
            
            document.addEventListener('mousemove', slideMove);
            document.addEventListener('mouseup', stopSlide);
        }
        
        function slideMove(e) {
            e.preventDefault();
            
            const pos = getCursorPos(e);
            if (pos < 0) pos = 0;
            if (pos > width) pos = width;
            
            slider.style.left = pos + 'px';
            beforeImg.style.width = pos + 'px';
        }
        
        function stopSlide() {
            document.removeEventListener('mousemove', slideMove);
            document.removeEventListener('mouseup', stopSlide);
        }
        
        function getCursorPos(e) {
            const rect = container.getBoundingClientRect();
            return e.pageX - rect.left - window.pageXOffset;
        }
    });
}

// Create and attach image comparison component
function createImageComparison(beforeSrc, afterSrc, container, label1, label2) {
    if (!container) return;
    
    const compContainer = document.createElement('div');
    compContainer.className = 'img-comp-container';
    
    // After image (background)
    const afterDiv = document.createElement('div');
    afterDiv.className = 'img-comp-img';
    const afterImg = document.createElement('img');
    afterImg.src = afterSrc;
    afterImg.width = "100%";
    const afterLabel = document.createElement('div');
    afterLabel.className = 'img-label img-after-label';
    afterLabel.textContent = label2 || 'After';
    afterDiv.appendChild(afterImg);
    afterDiv.appendChild(afterLabel);
    
    // Before image (foreground)
    const beforeDiv = document.createElement('div');
    beforeDiv.className = 'img-comp-img img-comp-before';
    const beforeImg = document.createElement('img');
    beforeImg.src = beforeSrc;
    beforeImg.width = "100%";
    const beforeLabel = document.createElement('div');
    beforeLabel.className = 'img-label img-before-label';
    beforeLabel.textContent = label1 || 'Before';
    beforeDiv.appendChild(beforeImg);
    beforeDiv.appendChild(beforeLabel);
    
    // Slider
    const slider = document.createElement('div');
    slider.className = 'img-comp-slider';
    
    // Assemble
    compContainer.appendChild(afterDiv);
    compContainer.appendChild(beforeDiv);
    compContainer.appendChild(slider);
    container.appendChild(compContainer);
    
    // Initialize
    setTimeout(function() {
        initComparisonSlider();
    }, 100);
}
