import React, { useEffect, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import Papa from "papaparse";
import "../styles/main.css";

const COLUMNS_TO_SHOW = ["cmpd_id", "SMILES"];

const CsvPreview = () => {
  const location = useLocation();
  const navigate = useNavigate();
  const baseUrl = "http://localhost:5001";
  const fileUrl = location.state?.fileUrl ? `${baseUrl}${location.state.fileUrl}` : null;

  const [rawData, setRawData] = useState([]);
  const [selectedRowIndex, setSelectedRowIndex] = useState(null);
  const [hoveredRowIndex, setHoveredRowIndex] = useState(null);

  useEffect(() => {
    if (fileUrl) {
      fetch(fileUrl)
        .then((response) => response.text())
        .then((csvText) => {
          const parsedResult = Papa.parse(csvText, { header: true });
          const filteredData = parsedResult.data.filter(
            (row) => Object.keys(row).length > 1
          );
          setRawData(filteredData);
        })
        .catch((error) => console.error("Error loading CSV:", error));
    }
  }, [fileUrl]);

  const handleRowClick = (row, index) => {
    const filename = location.state?.fileUrl?.split("/").pop();
    setSelectedRowIndex(index);
    navigate(`/details/${index}`, { state: { data: row,filename,totalRows: rawData.length, allData: rawData } });

  };

  return (
    <div className="app-container">
      {/* Header */}
      <header className="app-header">
        <div className="user-info" onClick={() => navigate("/")}> 
          <button className="avatar-button"> 
            <img src="/assets/home.png" alt="Home" className="home-icon" />
            <span className="home-label" style={{ color: "black" }}>Homepage</span>
          </button>
   
        </div>
      </header>

      {/* Main Section */}
      <main className="app-main">
        <h1 className="main-title" style={{ marginBottom: "100px" }}>CSV Preview</h1>
        {rawData.length > 0 ? (
          <div className="csv-container" style={{ maxHeight: "500px", maxWidth: "1000px", overflowY: "auto" }}>
            <table className="csv-table">
              <thead>
                <tr>
                  {COLUMNS_TO_SHOW.map((key) => (
                    <th key={key} style={{ color: "black" }}>{key}</th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {rawData.map((row, index) => (
                  <tr 
                    key={index} 
                    onClick={() => handleRowClick(row, index)}
                    onMouseEnter={() => setHoveredRowIndex(index)}
                    onMouseLeave={() => setHoveredRowIndex(null)}
                    style={{ 
                      cursor: "pointer", 
                      backgroundColor: selectedRowIndex === index ? "##ffaa00" : hoveredRowIndex === index ? "#ffaa00" : "white" 
                    }}
                  >
                    {COLUMNS_TO_SHOW.map((colKey) => (
                      <td key={colKey} style={{ color: "black" }}>{row[colKey]}</td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        ) : (
          <p style={{ color: "black" }}>Loading CSV...</p>
        )}
      </main>
    </div>
  );
};

export default CsvPreview;
