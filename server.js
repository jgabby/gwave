const express = require('express');
const app= express();
const gwave = require('./gwave.js');

app.get('/', (req, res) => {
    // console.log(req);
    res.writeHead(200, {"Content-Type": "application/json"});
    let result;
    try {
        console.log(req.query);
        let frequency = Number(req.query.frequency);
        result = {result: "success", chart: gwave.gw_chart(frequency)};
    } catch(e) {
        console.error(e);
        result = {result: "error", error: e};
    }
    res.end(JSON.stringify(result));
});

app.listen(3000,() => {
    console.log("Server started on Port 3000");
});